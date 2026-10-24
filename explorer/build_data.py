"""Build the explorer's data from benchmark results.

Writes `data.js` (complexes, ligand properties and every arm's per-complex outcome for one
engine) and one `poses/<engine>/<id>.js` per complex (crystal ligand, protein,
binding pocket and the docked poses of every arm), which the page loads when a complex is opened.

Usage:
    python explorer/build_data.py [--results DIR] [--engine vina] [--out explorer]
                                  [--share-url URL]

Runs are read from `<results>/<engine>/<arm>/results.csv` and `poses/`, as written by
`dockprepbench run`. Seed replicates (`baseline_seed<N>`) are not shown.
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
from dockprepbench.analysis.stats import load_run  # noqa: E402
from dockprepbench.datasets.posebusters import SET_NAME, download, load_complexes  # noqa: E402

ADDITIVES = {
    "HOH", "EDO", "GOL", "SO4", "PO4", "DMS", "PEG", "PG4", "ACT", "CL", "NA", "K", "BR",
    "IOD", "FMT", "MPD", "TRS", "EPE", "MES", "BME", "IMD", "NO3", "SCN", "1PE",
}  # fmt: skip
ARMS = {
    "baseline": ("Baseline", "4 Å padding, supplied ligand"),
    "box_padding_2": ("Box 2 Å", "Box padding 2 Å"),
    "box_padding_6": ("Box 6 Å", "Box padding 6 Å"),
    "box_padding_8": ("Box 8 Å", "Box padding 8 Å"),
    "ligand_dimorphite": ("Protonated ligand", "Ligand protonated at pH 7.4"),
    "ligand_rdkit_conformer": ("RDKit conformer", "Start conformer regenerated"),
}
ENGINES = {
    "vina": "AutoDock Vina 1.2.7",
    "smina": "Smina",
    "gnina": "Gnina 1.3.3",
    "rdock": "rDock 24.03",
}
POCKET_CUTOFF = 8.0  # Å from any crystal-ligand atom


def depiction(mol: Chem.Mol) -> str:
    """
    Draw a ligand in 2D as inline SVG, in `currentColor` so it follows the page theme.

    Args:
        mol (Chem.Mol): Ligand.

    Returns:
        str: Compact SVG markup.
    """
    flat = Chem.Mol(mol)
    rdDepictor.Compute2DCoords(flat)
    drawer = rdMolDraw2D.MolDraw2DSVG(240, 170, -1, -1, True)  # noFreetype: <text>, not glyph paths
    opts = drawer.drawOptions()
    opts.clearBackground = False
    opts.useBWAtomPalette()
    opts.padding = 0.06
    drawer.DrawMolecule(flat)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()
    svg = svg[svg.index("<svg") :]
    for colour in ("#000000", "#FF0000", "#000"):
        svg = svg.replace(colour, "currentColor")
    return compact_svg(svg)


def compact_svg(svg: str) -> str:
    """
    Move RDKit's repeated per-element styles onto the root to shrink each drawing ~4x.

    Args:
        svg (str): SVG from RDKit's drawer.

    Returns:
        str: Equivalent, smaller SVG.
    """

    def path(m):
        style = m.group(2)
        filled = "fill:currentColor" in style
        attrs = re.sub(r" class='[^']*'| style='[^']*'", "", m.group(1) + m.group(3))
        fill = " fill='currentColor'" if filled else ""
        return f"<path{attrs}{fill}"

    svg = re.sub(r"<path([^>]*?) style='([^']*)'([^>]*?)", path, svg)
    svg = re.sub(r"(<path[^>]*?) class='[^']*'", r"\1", svg)

    def text(m):
        anchor = re.search(r"text-anchor:(\w+)", m.group(0))
        size = re.search(r"font-size:([\d.]+)px", m.group(0))
        align = f" text-anchor='{anchor.group(1)}'" if anchor and anchor.group(1) != "start" else ""
        font_size = size.group(1) if size else 15
        return f"<text x='{m.group(1)}' y='{m.group(2)}' font-size='{font_size}'{align}>"

    svg = re.sub(r"<text x='([^']*)' y='([^']*)'[^>]*>", text, svg)
    svg = re.sub(r"<!--.*?-->", "", svg)
    svg = re.sub(r"<rect[^>]*/>", "", svg)
    svg = re.sub(
        r"<svg[^>]*?viewBox='([^']*)'[^>]*>",
        r"<svg xmlns='http://www.w3.org/2000/svg' viewBox='\1' fill='none' stroke='currentColor' "
        r"stroke-width='1.6' stroke-linecap='round' font-family='sans-serif'>",
        svg,
        count=1,
    )
    svg = svg.replace("<text ", "<text fill='currentColor' stroke='none' ")
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


def main() -> None:
    """Build data.js and the per-complex pose files."""
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--results", type=Path, default=None, help="Results root")
    parser.add_argument("--engine", default="vina", help="Engine to show")
    parser.add_argument("--out", type=Path, default=REPO / "explorer", help="Explorer directory")
    parser.add_argument("--share-url", default=None, help="Base URL for 'Copy link'")
    args = parser.parse_args()

    results = args.results or default_results()
    engine_dir = results / args.engine
    if not engine_dir.is_dir():
        sys.exit(f"No results for engine '{args.engine}' under {results}")
    arms = find_arms(engine_dir)
    if "baseline" not in arms:
        sys.exit(f"No finished baseline under {engine_dir}; nothing to compare with")
    runs = {arm: load_run(engine_dir / arm / "results.csv") for arm in arms}

    pose_dir = args.out / "poses" / args.engine
    shutil.rmtree(pose_dir, ignore_errors=True)
    pose_dir.mkdir(parents=True)

    dataset = REPO / "data/raw/posebusters"
    if not (dataset / SET_NAME).is_dir():
        print("Downloading the PoseBusters set (about 50 MB, once)...")
        dataset = download(REPO / "data/raw")

    rows = []
    for c in load_complexes(dataset):
        if c.id not in runs["baseline"].index:
            continue
        lig = next(Chem.SDMolSupplier(str(c.crystal_ligand)))
        hets = {
            ln[17:20].strip()
            for ln in c.protein.read_text().splitlines()
            if ln.startswith("HETATM")
        }
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
                "arms": {arm: outcome(run.loc[c.id]) for arm, run in runs.items()},
            }
        )
        payload = {
            "crystal": Chem.MolToMolBlock(lig),
            "pocket": pocket(c.protein, lig),
            "protein": chains(c.protein),
            "poses": {arm: poses(engine_dir / arm / "poses" / f"{c.id}.sdf") for arm in arms},
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
        "dataset": "PoseBusters Benchmark set",
        "engine": ENGINES.get(args.engine, args.engine),
        "arms": [
            {"key": a, "label": ARMS.get(a, (a, a))[0], "detail": ARMS.get(a, (a, a))[1]}
            for a in arms
        ],
        "complexes": rows,
        "real": {},
        "poseDir": f"poses/{args.engine}",
    }
    (args.out / "data.js").write_text(
        "window.EXPLORER_DATA = " + json.dumps(data, separators=(",", ":")) + ";\n",
        encoding="utf-8",
    )
    size = sum(f.stat().st_size for f in pose_dir.iterdir()) / 2**20
    print(
        f"{len(rows)} complexes, {len(arms)} arms ({', '.join(arms)}) for {args.engine}; "
        f"data.js {(args.out / 'data.js').stat().st_size / 2**10:.0f} KB, poses {size:.0f} MB"
    )


if __name__ == "__main__":
    main()
