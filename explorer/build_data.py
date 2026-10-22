"""Build data.js for the explorer prototype.

Real: complex IDs, ligand properties, cofactor flags, and poses/pockets for complexes that
have a smoke-test run. Demo: per-arm outcomes (seeded, flagged `demo: true`).

Usage: python explorer/build_data.py REPO_ROOT OUT_DIR [SHARE_URL]

SHARE_URL is the base of the "Copy link" button (e.g. the hosted page); without it the page
uses its own address.
"""

import json
import re
import sys
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors, rdDepictor, rdMolDescriptors
from rdkit.Chem.Draw import rdMolDraw2D

repo, out_dir = Path(sys.argv[1]), Path(sys.argv[2])
share_url = sys.argv[3] if len(sys.argv) > 3 else None
sys.path.insert(0, str(repo / "src"))
from dockprepbench.datasets.posebusters import load_complexes  # noqa: E402

ADDITIVES = {
    "HOH",
    "EDO",
    "GOL",
    "SO4",
    "PO4",
    "DMS",
    "PEG",
    "PG4",
    "ACT",
    "CL",
    "NA",
    "K",
    "BR",
    "IOD",
    "FMT",
    "MPD",
    "TRS",
    "EPE",
    "MES",
    "BME",
    "IMD",
    "NO3",
    "SCN",
    "1PE",
}
ARMS = [
    ("baseline", "Baseline", "4 Å padding, supplied ligand"),
    ("box_padding_2", "Box 2 Å", "Box padding 2 Å"),
    ("box_padding_6", "Box 6 Å", "Box padding 6 Å"),
    ("box_padding_8", "Box 8 Å", "Box padding 8 Å"),
    ("ligand_dimorphite", "Protonated ligand", "Ligand protonated at pH 7.4"),
    ("ligand_rdkit_conformer", "RDKit conformer", "Start conformer regenerated"),
]


def depiction(mol: Chem.Mol) -> str:
    """2D drawing as inline SVG, drawn in `currentColor` so it follows the page theme."""
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
    """Move RDKit's repeated per-element styles onto the root to shrink each drawing ~4x."""

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


complexes = load_complexes(repo / "data/raw/posebusters")
rows = []
for c in complexes:
    lig = next(Chem.SDMolSupplier(str(c.crystal_ligand)))
    hets = {
        line[17:20].strip()
        for line in c.protein.read_text().splitlines()
        if line.startswith("HETATM")
    }
    cofactors = sorted(hets - ADDITIVES)
    rows.append(
        {
            "id": c.id,
            "heavy": lig.GetNumHeavyAtoms(),
            "rotb": rdMolDescriptors.CalcNumRotatableBonds(lig),
            "mw": round(Descriptors.MolWt(lig), 1),
            "cofactors": cofactors[:4],
            "svg": depiction(lig),
        }
    )

# Demo outcomes: success odds fall with flexibility; arms nudge them; seeded and paired.
rng = np.random.default_rng(7)
latent = rng.normal(0, 1, len(rows))
arm_shift = {
    "baseline": 0.0,
    "box_padding_2": 0.25,
    "box_padding_6": -0.15,
    "box_padding_8": -0.4,
    "ligand_dimorphite": 0.1,
    "ligand_rdkit_conformer": -0.05,
}
for row, z in zip(rows, latent, strict=True):
    row["arms"] = {}
    for key, _, _ in ARMS:
        score = 0.3 + z - 0.08 * row["rotb"] + arm_shift[key] + rng.normal(0, 0.45)
        ok = score > 0
        top1 = float(rng.gamma(2, 0.4)) if ok else float(rng.uniform(2.2, 11))
        best = (
            top1
            if ok
            else (float(rng.gamma(2, 0.45)) if rng.random() < 0.35 else top1 * rng.uniform(0.5, 1))
        )
        cof_clash = bool(row["cofactors"]) and ok and rng.random() < 0.15
        row["arms"][key] = {
            "rmsd": round(top1, 2),
            "best": round(min(best, top1), 2),
            "valid": bool(ok and not cof_clash and rng.random() > 0.08),
            "shift": round(top1 * rng.uniform(0.3, 0.9), 2),
            "failed": ["volume_overlap_with_organic_cofactors"] if cof_clash else [],
        }

# Real baseline outcomes override the demo ones wherever a real run exists.
import pandas as pd  # noqa: E402

smoke = pd.read_csv(repo / "results/smoke_poses/results.csv").set_index("id")
by_id = {row["id"]: row for row in rows}
for cid, r in smoke.iterrows():
    by_id[cid]["arms"]["baseline"] = {
        "rmsd": round(float(r["top1_rmsd"]), 2),
        "best": round(float(r["best_rmsd"]), 2),
        "valid": bool(r["top1_pb_valid"]),
        "shift": round(float(r["top1_centroid_shift"]), 2),
        "failed": [c for c in str(r["top1_failed_checks"]).split(";") if c and c != "nan"],
        "real": True,
    }

# Real poses: any complex with a saved smoke-test run.
real = {}
pose_dir = repo / "results/smoke_poses/poses"
for sdf in sorted(pose_dir.glob("*.sdf")):
    cid = sdf.stem
    c = next(x for x in complexes if x.id == cid)
    crystal = next(Chem.SDMolSupplier(str(c.crystal_ligand)))
    centre = crystal.GetConformer().GetPositions()
    receptor = (repo / f"data/prepared/ph7.4/{cid}_receptor.pdb").read_text().splitlines()
    near = set()
    for line in receptor:
        if line.startswith("ATOM"):
            xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            if np.min(np.linalg.norm(centre - xyz, axis=1)) < 8.0:
                near.add(line[17:27])
    pocket = "\n".join(ln for ln in receptor if ln.startswith("ATOM") and ln[17:27] in near)
    poses = [m for m in Chem.SDMolSupplier(str(sdf), removeHs=False)]
    real[cid] = {
        "crystal": Chem.MolToMolBlock(crystal),
        "pocket": pocket + "\nEND\n",
        "poses": [
            {
                "rank": m.GetIntProp("rank"),
                "score": round(m.GetDoubleProp("score"), 2),
                "rmsd": round(m.GetDoubleProp("rmsd"), 2),
                "block": Chem.MolToMolBlock(Chem.RemoveHs(m)),
            }
            for m in poses
        ],
    }

data = {
    "demo": True,
    "shareUrl": share_url,
    "dataset": "PoseBusters Benchmark set",
    "engine": "AutoDock Vina 1.2.7",
    "arms": [{"key": key, "label": label, "detail": detail} for key, label, detail in ARMS],
    "complexes": rows,
    "real": real,
}
out_dir.mkdir(parents=True, exist_ok=True)
(out_dir / "data.js").write_text(
    "window.EXPLORER_DATA = " + json.dumps(data, separators=(",", ":")) + ";\n", encoding="utf-8"
)
print(
    len(rows),
    "complexes;",
    "real poses for",
    list(real),
    ";",
    round((out_dir / "data.js").stat().st_size / 1024),
    "KB",
)
