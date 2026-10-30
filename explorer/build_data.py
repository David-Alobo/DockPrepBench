"""Build the explorer's data from benchmark results.

Writes `data.js` (complexes, ligand properties and the per-complex outcome of every engine and
arm) and one `poses/<id>.js` per complex (crystal ligand, protein, binding pocket and the docked
poses of every engine and arm), which the page loads when a complex is opened.

Usage:
    python explorer/build_data.py [--results DIR] [--engine NAME ...] [--out explorer]
                                  [--share-url URL]

Runs are read from `<results>/<engine>/<arm>/results.csv` and `poses/`, as written by
`dockprepbench run`. Every engine with a finished baseline is included unless `--engine` names
some; each shows the arms it has finished. Seed replicates (`baseline_seed<N>`) are not shown.
"""

import argparse
import json
import re
import shutil
import sys
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors, rdDepictor, rdMolDescriptors
from rdkit.Chem.Draw import rdMolDraw2D

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
from dockprepbench.analysis.agreement import collect as collect_agreement  # noqa: E402
from dockprepbench.analysis.sensitivity import ADDITIVES  # noqa: E402
from dockprepbench.analysis.stats import load_run  # noqa: E402
from dockprepbench.datasets.posebusters import SET_NAME, SETS, download, load_complexes  # noqa: E402

# Dataset key -> (label shown on the page, folder in the PoseBusters archive, file suffix).
DATASETS = {
    "posebusters": ("PoseBusters Benchmark set", SET_NAME, ""),
    "astex": ("Astex Diverse Set", SETS["astex"], "_astex"),
}

ARMS = {
    "baseline": ("Baseline", "4 Å padding, supplied ligand"),
    "box_padding_2": ("Box 2 Å", "Box padding 2 Å"),
    "box_padding_6": ("Box 6 Å", "Box padding 6 Å"),
    "box_padding_8": ("Box 8 Å", "Box padding 8 Å"),
    "ligand_dimorphite": ("Protonated ligand", "Ligand protonated at pH 7.4"),
    "ligand_rdkit_conformer": ("RDKit conformer", "Start conformer regenerated"),
    # Follow-up arms (Vina only): protonation combined with each padding, and a box sized to the
    # ligand.
    "ligand_dimorphite_box_2": ("Protonated + 2 Å", "Ligand protonated at pH 7.4, box padding 2 Å"),
    "ligand_dimorphite_box_6": ("Protonated + 6 Å", "Ligand protonated at pH 7.4, box padding 6 Å"),
    "ligand_dimorphite_box_8": ("Protonated + 8 Å", "Ligand protonated at pH 7.4, box padding 8 Å"),
    "box_gyration": ("Box from Rg", "Cube 2.857 times the ligand's radius of gyration"),
}
ENGINES = {
    "vina": "AutoDock Vina 1.2.7",
    "smina": "Smina",
    "gnina": "Gnina 1.3.3",
    "rdock": "rDock 24.03",
}
POCKET_CUTOFF = 5.0  # Å from any crystal-ligand atom: residues in contact with it


DEPICTION_SIZE = (300, 120)  # px; wide, like the header box it is shown in
# RDKit's default atom colours -> colours legible on light and dark backgrounds. Carbon (black)
# becomes currentColor, so the drawing follows the page's text colour.
ELEMENT_COLOURS = {
    "#000000": "currentColor",
    "#0000FF": "#3b6fe0",  # N
    "#FF0000": "#e03131",  # O
    "#CCCC00": "#d4a800",  # S
    "#33CCCC": "#1fa2a2",  # F
    "#00CC00": "#2f9e44",  # Cl
    "#FF7F00": "#e8590c",  # P
    "#A62929": "#b5542a",  # Br
    "#940094": "#9a4fd1",  # I
}


def depiction(mol: Chem.Mol) -> str:
    """
    Draw a ligand in 2D as inline SVG, carbon in `currentColor` and heteroatoms coloured.

    Args:
        mol (Chem.Mol): Ligand.

    Returns:
        str: Compact SVG markup.
    """
    flat = Chem.Mol(mol)
    rdDepictor.Compute2DCoords(flat)
    drawer = rdMolDraw2D.MolDraw2DSVG(*DEPICTION_SIZE, -1, -1, True)  # noFreetype: <text> labels
    opts = drawer.drawOptions()
    opts.clearBackground = False
    opts.padding = 0.04
    opts.bondLineWidth = 2.5
    drawer.DrawMolecule(flat)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()
    svg = svg[svg.index("<svg") :]
    svg = re.sub(
        r"#[0-9A-Fa-f]{6}",
        lambda m: ELEMENT_COLOURS.get(m.group(0).upper(), m.group(0)),
        svg,
    )
    return compact_svg(svg)


def compact_svg(svg: str) -> str:
    """
    Move RDKit's repeated per-element styles onto the root to shrink each drawing ~4x.

    Keeps each element's stroke and fill colour when it differs from the root's.

    Args:
        svg (str): SVG from RDKit's drawer.

    Returns:
        str: Equivalent, smaller SVG.
    """

    def colour(style: str, prop: str) -> str | None:
        found = re.search(rf"(?:^|;){prop}:([^;]+)", style)
        return found.group(1) if found else None

    def path(m):
        style = m.group(2)
        attrs = re.sub(r" class='[^']*'| style='[^']*'", "", m.group(1) + m.group(3))
        stroke, fill = colour(style, "stroke"), colour(style, "fill")
        if stroke and stroke not in ("currentColor", "none"):
            attrs += f" stroke='{stroke}'"
        if fill and fill != "none":
            attrs += f" fill='{fill}'"
        return f"<path{attrs}"

    svg = re.sub(r"<path([^>]*?) style='([^']*)'([^>]*?)", path, svg)
    svg = re.sub(r"(<path[^>]*?) class='[^']*'", r"\1", svg)

    def text(m):
        style = m.group(0)
        anchor = re.search(r"text-anchor:(\w+)", style)
        size = re.search(r"font-size:([\d.]+)px", style)
        align = f" text-anchor='{anchor.group(1)}'" if anchor and anchor.group(1) != "start" else ""
        font_size = size.group(1) if size else 15
        fill = colour(re.search(r"style='([^']*)'", style).group(1), "fill") or "currentColor"
        return (
            f"<text x='{m.group(1)}' y='{m.group(2)}' font-size='{font_size}' fill='{fill}' "
            f"stroke='none'{align}>"
        )

    svg = re.sub(r"<text x='([^']*)' y='([^']*)'[^>]*>", text, svg)
    svg = re.sub(r"<!--.*?-->", "", svg)
    svg = re.sub(r"<rect[^>]*/>", "", svg)
    svg = re.sub(
        r"<svg[^>]*?viewBox='([^']*)'[^>]*>",
        r"<svg xmlns='http://www.w3.org/2000/svg' viewBox='\1' fill='none' stroke='currentColor' "
        r"stroke-width='2.5' stroke-linecap='round' font-family='sans-serif'>",
        svg,
        count=1,
    )
    return " ".join(svg.split())


def default_results() -> Path:
    """
    Pick the results directory: the local copy of the server results if present.

    Returns:
        Path: `results/server/results` when it exists, otherwise `results`.
    """
    server = REPO / "results/server/results"
    return server if server.is_dir() else REPO / "results"


def find_arms(engine_dir: Path) -> list[str]:
    """
    List an engine's finished arms, known arms first in their usual order.

    Args:
        engine_dir (Path): `<results>/<engine>`.

    Returns:
        list[str]: Arm names with a `results.csv`, excluding seed replicates.
    """
    done = {
        d.name
        for d in engine_dir.iterdir()
        if (d / "results.csv").is_file() and not re.match(r"baseline_seed\d+$", d.name)
    }
    return [a for a in ARMS if a in done] + sorted(done - set(ARMS))


def outcome(row) -> dict:
    """
    Turn one results row into the page's outcome record.

    Args:
        row (pd.Series): Row of a run loaded with `load_run`.

    Returns:
        dict: rmsd, best, valid, shift and failed checks; or `error` if the complex failed.
    """
    if row["error"]:
        return {"error": row["error"].splitlines()[0][:200]}
    # RMSDs keep full precision: rounding 1.996 Å to 2.00 would turn a success into a miss.
    return {
        "rmsd": float(row["top1_rmsd"]),
        "best": float(row["best_rmsd"]),
        "valid": bool(row["top1_pb_valid"]),
        "shift": round(float(row["top1_centroid_shift"]), 2),
        "failed": [c for c in row["top1_failed_checks"].split(";") if c],
    }


def pocket(protein: Path, ligand: Chem.Mol) -> str:
    """
    Cut the binding pocket out of the dataset's protein file, cofactors included.

    Args:
        protein (Path): Protein PDB as distributed with the dataset.
        ligand (Chem.Mol): Crystal ligand, which defines the site.

    Returns:
        str: PDB records of every residue with an atom within POCKET_CUTOFF of the ligand.
    """
    atoms = [
        ln
        for ln in protein.read_text().splitlines()
        if ln.startswith(("ATOM", "HETATM")) and ln[17:20].strip() != "HOH"
    ]
    xyz = np.array([[float(ln[30:38]), float(ln[38:46]), float(ln[46:54])] for ln in atoms])
    site = ligand.GetConformer().GetPositions()
    dist = np.linalg.norm(xyz[:, None, :] - site[None, :, :], axis=2).min(axis=1)
    near = {ln[17:27] for ln, d in zip(atoms, dist, strict=True) if d < POCKET_CUTOFF}
    return "\n".join(ln for ln in atoms if ln[17:27] in near) + "\nEND\n"


def chains(protein: Path) -> str:
    """
    Keep the protein's heavy atoms, from which the page draws a cartoon and a surface.

    Args:
        protein (Path): Protein PDB as distributed with the dataset.

    Returns:
        str: Heavy-atom ATOM records, trimmed after the coordinates to save space.
    """
    lines = [
        ln[:54]
        for ln in protein.read_text().splitlines()
        if ln.startswith("ATOM")
        and ln[76:78].strip() != "H"
        and not ln[12:16].strip().startswith("H")
    ]
    return "\n".join(lines) + "\nEND\n"


def poses(sdf: Path) -> list[dict]:
    """
    Read one complex's saved poses.

    Args:
        sdf (Path): `<run>/poses/<id>.sdf`.

    Returns:
        list[dict]: rank, score, RMSD and heavy-atom mol block per pose; empty if missing.
    """
    if not sdf.is_file():
        return []
    return [
        {
            "rank": m.GetIntProp("rank"),
            "score": round(m.GetDoubleProp("score"), 2),
            "rmsd": m.GetDoubleProp("rmsd"),
            "block": Chem.MolToMolBlock(Chem.RemoveHs(m, sanitize=False)),
        }
        for m in Chem.SDMolSupplier(str(sdf), removeHs=False, sanitize=False)
        if m is not None
    ]


def find_engines(results: Path, wanted: list[str] | None) -> list[str]:
    """
    List the engines to show: those requested, or every engine with a finished baseline.

    Args:
        results (Path): Results root.
        wanted (list[str] | None): Engines named on the command line, or None for all.

    Returns:
        list[str]: Engine names, known engines first in their usual order.
    """
    if wanted:
        missing = [e for e in wanted if not (results / e / "baseline" / "results.csv").is_file()]
        if missing:
            sys.exit(f"No finished baseline under {results} for: {', '.join(missing)}")
        return wanted
    found = {d.name for d in results.iterdir() if (d / "baseline" / "results.csv").is_file()}
    if not found:
        sys.exit(f"No engine under {results} has a finished baseline; nothing to show")
    return [e for e in ENGINES if e in found] + sorted(found - set(ENGINES))


def main() -> None:
    """Build data.js and the per-complex pose files."""
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--results", type=Path, default=None, help="Results root")
    parser.add_argument(
        "--engine", action="append", help="Engine to include (repeatable); default: all"
    )
    parser.add_argument("--out", type=Path, default=REPO / "explorer", help="Explorer directory")
    parser.add_argument("--share-url", default=None, help="Base URL for 'Copy link'")
    parser.add_argument(
        "--dataset",
        choices=sorted(DATASETS),
        default="posebusters",
        help="Complex set; astex writes data_astex.js and poses_astex/ next to the default",
    )
    args = parser.parse_args()
    label, set_name, suffix = DATASETS[args.dataset]

    results = args.results or (
        default_results() if args.dataset == "posebusters" else default_results().parent / "astex"
    )
    engines = find_engines(results, args.engine)
    arms = {e: find_arms(results / e) for e in engines}
    runs = {
        e: {arm: load_run(results / e / arm / "results.csv") for arm in arms[e]} for e in engines
    }

    # Agreement of each baseline top pose with the primary arms (paper section on agreement).
    print("Computing pose agreement across arms...")
    agreement = collect_agreement(results, engines)
    agree = {(r.engine, r.id): int(r.prep_agree) for r in agreement.itertuples()}

    pose_dir = args.out / f"poses{suffix}"
    shutil.rmtree(pose_dir, ignore_errors=True)
    pose_dir.mkdir(parents=True)

    dataset = REPO / "data/raw/posebusters"
    if not (dataset / SET_NAME).is_dir():
        print("Downloading the PoseBusters set (about 50 MB, once)...")
        dataset = download(REPO / "data/raw")

    rows = []
    for c in load_complexes(dataset, set_name):
        lig = next(Chem.SDMolSupplier(str(c.crystal_ligand)))
        records = c.protein.read_text().splitlines()
        hets = {ln[17:20].strip() for ln in records if ln.startswith("HETATM")}
        residues = {ln[21:27] for ln in records if ln.startswith("ATOM")}  # chain, number, icode
        # Other copies of the ligand of interest (its CCD code ends the ID) are not cofactors.
        cofactors = sorted(hets - ADDITIVES - {c.id.split("_")[-1]})
        rows.append(
            {
                "id": c.id,
                "heavy": lig.GetNumHeavyAtoms(),
                "rotb": rdMolDescriptors.CalcNumRotatableBonds(lig),
                "mw": round(Descriptors.MolWt(lig), 1),
                "cofactors": cofactors[:4],
                "svg": depiction(lig),
                "smiles": Chem.MolToSmiles(Chem.RemoveHs(lig)),
                "charge": Chem.GetFormalCharge(lig),
                "chains": len({r[0] for r in residues}),
                "residues": len(residues),
                # A complex missing from a run (e.g. a run limited with --limit) counts as failed.
                "byEngine": {
                    e: {
                        arm: outcome(run.loc[c.id]) if c.id in run.index else {"error": "not run"}
                        for arm, run in runs[e].items()
                    }
                    for e in engines
                },
                "agreeByEngine": {e: agree[(e, c.id)] for e in engines if (e, c.id) in agree},
            }
        )
        payload = {
            "crystal": Chem.MolToMolBlock(lig),
            "pocket": pocket(c.protein, lig),
            "protein": chains(c.protein),
            "poses": {
                e: {arm: poses(results / e / arm / "poses" / f"{c.id}.sdf") for arm in arms[e]}
                for e in engines
            },
        }
        (pose_dir / f"{c.id}.js").write_text(
            f"window.ExplorerPoses({json.dumps(c.id)}, "
            + json.dumps(payload, separators=(",", ":"))
            + ");\n",
            encoding="utf-8",
        )

    data = {
        "demo": False,
        "shareUrl": args.share_url,
        "dataset": label,
        "datasetKey": args.dataset,
        "engines": [
            {
                "key": e,
                "label": ENGINES.get(e, e),
                "arms": [
                    {"key": a, "label": ARMS.get(a, (a, a))[0], "detail": ARMS.get(a, (a, a))[1]}
                    for a in arms[e]
                ],
            }
            for e in engines
        ],
        "complexes": rows,
        "real": {},
        "poseDir": f"poses{suffix}",
    }
    (args.out / f"data{suffix}.js").write_text(
        "window.EXPLORER_DATA = " + json.dumps(data, separators=(",", ":")) + ";\n",
        encoding="utf-8",
    )
    size = sum(f.stat().st_size for f in pose_dir.iterdir()) / 2**20
    print(
        f"{len(rows)} complexes; data{suffix}.js "
        f"{(args.out / f'data{suffix}.js').stat().st_size / 2**10:.0f} KB, "
        f"poses {size:.0f} MB"
    )
    for e in engines:
        print(f"  {e}: {len(arms[e])} arms ({', '.join(arms[e])})")


if __name__ == "__main__":
    main()
