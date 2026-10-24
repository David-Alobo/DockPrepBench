"""Draw the study-design figure (Figure 1) as SVG.

Usage:
    python scripts/make_concept_figure.py [OUT.svg]                       # version 2
    python scripts/make_concept_figure.py --renders docs/figures/renders  # version 3
    python scripts/make_concept_figure.py --renders docs/figures/renders_v4 --v4  # version 4

Version 2 (docs/figures/concept_v2.svg) shows the study as run: the receptor is held fixed,
and the arms vary the search box and the ligand; panel a uses drawn glyphs. Version 3
(docs/figures/concept_v3.svg) is the same figure with panel a showing a real complex: 3D
renders made by scripts/figure_renders/render.ps1 and the ligand's 2D structure drawn by
RDKit. Version 1 (concept.svg, the original plan with receptor, cofactor and water arms) is
kept for reference.

The figure is 180 mm wide (double column). At that size 1 px = 0.32 pt, so body text
(22 px) prints at 7 pt. Convert to PDF/PNG with e.g. `rsvg-convert`.
"""

import argparse
import base64
import io
import math
from pathlib import Path

parser = argparse.ArgumentParser(description="Draw Figure 1 (study design) as SVG.")
parser.add_argument("out", nargs="?", type=Path, help="Output SVG (default by version)")
parser.add_argument("--renders", type=Path, help="Directory of 3D renders: makes version 3")
parser.add_argument("--complex", default="5SAK_ZRY", help="Complex shown in versions 3 and 4")
parser.add_argument(
    "--v4",
    action="store_true",
    help="Version 4 (with --renders from render.ps1 -Style v4): the ligand card shows the "
    "supplied and pH 7.4 states, and the evaluation picture is labelled with its RMSD",
)
ARGS = parser.parse_args()

WIDTH, HEIGHT = 1600, 1040
FONT = "Arial, Helvetica, sans-serif"

# Light-mode palette: ink for text, categorical slots for the stages (same as analysis/figures).
INK, INK_2, MUTED, GRID = "#0b0b0b", "#52514e", "#898781", "#e1e0d9"
SURFACE = "#ffffff"
RECEPTOR, LIGAND, BOX, ENGINE = "#2a78d6", "#eb6834", "#1baf7a", "#6250d6"
TINT = {RECEPTOR: "#eef4fc", LIGAND: "#fdf0ea", BOX: "#e9f7f1", ENGINE: "#efedfb", INK_2: "#f4f4f2"}

BODY, TITLE, PANEL = 22, 26, 32

out: list[str] = []


def add(element: str) -> None:
    out.append(element)


def text(x, y, s, size=BODY, color=INK_2, weight="normal", anchor="start", style="normal"):
    add(
        f'<text x="{x:.1f}" y="{y:.1f}" font-family="{FONT}" font-size="{size}" '
        f'font-weight="{weight}" font-style="{style}" fill="{color}" '
        f'text-anchor="{anchor}">{s}</text>'
    )


def hexagon(cx, cy, r, stroke, fill="none", width=2.5):
    pts = " ".join(
        f"{cx + r * math.cos(math.radians(a)):.1f},{cy + r * math.sin(math.radians(a)):.1f}"
        for a in range(0, 360, 60)
    )
    add(f'<polygon points="{pts}" fill="{fill}" stroke="{stroke}" stroke-width="{width}"/>')


def ligand(cx, cy, color, fill="none", r=13, dashed=False):
    """A ring with a short side chain: the generic small-molecule glyph."""
    dash = ' stroke-dasharray="5 4"' if dashed else ""
    pts = " ".join(
        f"{cx + r * math.cos(math.radians(a)):.1f},{cy + r * math.sin(math.radians(a)):.1f}"
        for a in range(0, 360, 60)
    )
    add(f'<polygon points="{pts}" fill="{fill}" stroke="{color}" stroke-width="2.5"{dash}/>')
    add(
        f'<polyline points="{cx + r:.1f},{cy:.1f} {cx + r + 12:.1f},{cy - 8:.1f} '
        f'{cx + r + 24:.1f},{cy:.1f}" fill="none" stroke="{color}" stroke-width="2.5"{dash}/>'
    )


def protein(cx, cy, stroke, fill):
    """Rounded blob with a binding pocket notched into its top edge."""
    add(
        f'<path d="M{cx - 62},{cy + 8} C{cx - 62},{cy - 30} {cx - 34},{cy - 40} '
        f"{cx - 20},{cy - 38} "
        f"L{cx - 16},{cy - 14} Q{cx},{cy - 2} {cx + 16},{cy - 14} L{cx + 20},{cy - 38} "
        f"C{cx + 42},{cy - 40} {cx + 62},{cy - 24} {cx + 62},{cy + 6} "
        f"C{cx + 62},{cy + 36} {cx + 30},{cy + 42} {cx},{cy + 40} "
        f'C{cx - 36},{cy + 42} {cx - 62},{cy + 34} {cx - 62},{cy + 8} Z" '
        f'fill="{fill}" stroke="{stroke}" stroke-width="2.5"/>'
    )


def arrow(x1, y1, x2, y2, color=MUTED, width=2.5):
    add(
        f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="{color}" '
        f'stroke-width="{width}" marker-end="url(#head)"/>'
    )


# Picture window in a card (versions 3 and 4): below the title, above the text.
PIC_W, PIC_H, PIC_TOP = 196, 108, 60
LIGAND_W = 210  # version 4 ligand drawing: the card's full width, less a small margin
# Version 4 search-space window: the card's full width and the height between title and text.
BOX_W, BOX_H, BOX_TOP = 210, 112, 60


def png_data(
    path: Path,
    size: tuple[int, int] = (PIC_W, PIC_H),
    fit_to: str | None = None,
    margin: float = 1.06,
) -> str:
    """
    Crop a render to its drawing and return it as a PNG data URI.

    Renders have a white background (headless WebGL exports no alpha), so the crop is to the
    pixels that differ from white, padded and fitted to the window's aspect ratio.

    Args:
        path (Path): Render from scripts/figure_renders/render.ps1.
        size (tuple[int, int]): Window width and height in figure px.
        fit_to (str | None): Hex colour to crop to instead of all the drawing, so an object
            of that colour fills the window and anything around it runs off the edges.
        margin (float): Crop size relative to the drawing, so above 1 leaves room around it.

    Returns:
        str: data:image/png;base64,... URI.
    """
    from PIL import Image, ImageChops

    image = Image.open(path).convert("RGB")
    ink = ImageChops.difference(image, Image.new("RGB", image.size, fit_to or "white"))
    if fit_to:
        # Pixels close to the colour: the difference from it is small in every channel (a
        # grey-level difference lets shaded ribbons through).
        r, g, b = ink.split()
        largest = ImageChops.lighter(ImageChops.lighter(r, g), b)
        mask = largest.point(lambda v: 255 if v < 50 else 0)
    else:
        mask = ink.point(lambda v: 255 if v > 12 else 0)
    left, top, right, bottom = mask.getbbox()
    # Grow the crop to the window's aspect ratio, centred on the drawing, with a margin.
    pic_w, pic_h = size
    cx, cy, w, h = (left + right) / 2, (top + bottom) / 2, right - left, bottom - top
    w, h = max(w, h * pic_w / pic_h) * margin, max(h, w * pic_h / pic_w) * margin
    box = (int(cx - w / 2), int(cy - h / 2), int(cx + w / 2), int(cy + h / 2))
    # Paste only the part inside the render, so any area beyond its edge stays white
    # (Image.crop would fill it black).
    canvas = Image.new("RGB", (box[2] - box[0], box[3] - box[1]), "white")
    inside = (max(box[0], 0), max(box[1], 0), min(box[2], image.width), min(box[3], image.height))
    canvas.paste(image.crop(inside), (inside[0] - box[0], inside[1] - box[1]))
    canvas = canvas.resize((pic_w * 3, pic_h * 3), Image.LANCZOS)  # 3x for print
    buffer = io.BytesIO()
    white_to_alpha(canvas).save(buffer, format="PNG", optimize=True)
    return "data:image/png;base64," + base64.b64encode(buffer.getvalue()).decode()


def white_to_alpha(image):
    """
    Make a picture drawn on white transparent, keeping antialiased edges clean.

    Each pixel becomes the colour that, laid over white at the derived opacity, reproduces it
    ("colour to alpha"): pure white turns fully transparent and the soft edges of ribbons and
    sticks blend into any background instead of leaving a white halo.

    Args:
        image (PIL.Image.Image): RGB picture on a white background.

    Returns:
        PIL.Image.Image: RGBA picture.
    """
    import numpy as np
    from PIL import Image

    rgb = np.asarray(image.convert("RGB"), dtype=float) / 255
    alpha = (1 - rgb).max(axis=2, keepdims=True)
    colour = 1 - (1 - rgb) / np.maximum(alpha, 1e-6)
    rgba = np.concatenate([np.clip(colour, 0, 1), alpha], axis=2)
    return Image.fromarray((rgba * 255).round().astype("uint8"), "RGBA")


def ligand_svg_data(complex_id: str) -> str:
    """
    Draw a complex's crystal ligand in 2D with RDKit, as an SVG document.

    Args:
        complex_id (str): PoseBusters complex id, e.g. 5SAK_ZRY.

    Returns:
        str: SVG document, heteroatoms in RDKit's colours.
    """
    from rdkit import Chem
    from rdkit.Chem import rdDepictor
    from rdkit.Chem.Draw import rdMolDraw2D

    sdf = Path("data/raw/posebusters/posebusters_benchmark_set") / complex_id
    mol = Chem.RemoveHs(next(Chem.SDMolSupplier(str(sdf / f"{complex_id}_ligand.sdf"))))
    rdDepictor.Compute2DCoords(mol)
    drawer = rdMolDraw2D.MolDraw2DSVG(PIC_W * 2, PIC_H * 2)
    opts = drawer.drawOptions()
    opts.bondLineWidth = 2.5
    opts.padding = 0.08
    opts.clearBackground = False  # transparent, so the card's tint shows through
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    return drawer.GetDrawingText()


def ligand_states_svg_data(complex_id: str) -> str:
    """
    Draw a ligand's dominant state at pH 7.4, highlighting the atoms protonation changes.

    One molecule at full size reads better in a card than the two states side by side; the
    legend names the change and the card text names both states. Uses the same protonation
    code as the benchmark.

    Args:
        complex_id (str): PoseBusters complex id.

    Returns:
        tuple[str, int]: SVG document, and the change in formal charge.
    """
    from rdkit import Chem
    from rdkit.Chem import rdDepictor
    from rdkit.Chem.Draw import rdMolDraw2D

    from dockprepbench.io import load_sdf_file
    from dockprepbench.preparation.ligand import protonate_ligand

    folder = Path("data/raw/posebusters/posebusters_benchmark_set") / complex_id
    supplied = load_sdf_file(folder / f"{complex_id}_ligand_start_conf.sdf")[0]
    # RemoveAllHs: RemoveHs keeps hydrogens that define double-bond stereo, so the two
    # states could differ in atom count.
    states = [Chem.RemoveAllHs(supplied), Chem.RemoveAllHs(protonate_ligand(supplied, 7.4))]
    changed = [
        i
        for i, (a, b) in enumerate(zip(states[0].GetAtoms(), states[1].GetAtoms(), strict=True))
        if a.GetFormalCharge() != b.GetFormalCharge()
    ]
    protonated = states[1]
    rdDepictor.Compute2DCoords(protonated)
    charge = Chem.GetFormalCharge(protonated) - Chem.GetFormalCharge(states[0])
    # Drawn at the size it is shown (1 unit = 1 figure px), so label sizes are print sizes:
    # 15 px is 4.8 pt at 180 mm, small enough to keep the atom labels clear of the bonds.
    # The charge label is set in the figure's own font, outside.
    drawer = rdMolDraw2D.MolDraw2DSVG(LIGAND_W, PIC_H)
    opts = drawer.drawOptions()
    opts.bondLineWidth = 1.6
    opts.padding = 0.03
    opts.minFontSize = 15
    opts.maxFontSize = 15
    opts.clearBackground = False
    drawer.DrawMolecule(
        protonated,
        highlightAtoms=changed,
        highlightAtomColors={i: (0.98, 0.72, 0.58) for i in changed},  # light ligand orange
        highlightBonds=[],
    )
    drawer.FinishDrawing()
    return drawer.GetDrawingText(), charge


def top_pose_rmsd(complex_id: str) -> float:
    """
    Read the RMSD of Vina's top baseline pose from the explorer's pose file for a complex.

    Args:
        complex_id (str): PoseBusters complex id.

    Returns:
        float: Heavy-atom RMSD in Å.
    """
    import json

    text = Path(f"explorer/poses/{complex_id}.js").read_text(encoding="utf-8")
    payload = json.loads(text[text.index(",") + 1 : text.rindex(")")])
    return next(p["rmsd"] for p in payload["poses"]["vina"]["baseline"] if p["rank"] == 1)


def inline_svg(svg: str, x: float, y: float, w: float, h: float) -> None:
    """
    Place an SVG drawing (e.g. from RDKit) in the figure as vector graphics.

    Embedding it as an <image> would make PDF converters rasterise it at low resolution.

    Args:
        svg (str): Complete SVG document.
        x (float): Left edge in figure units.
        y (float): Top edge in figure units.
        w (float): Width in figure units.
        h (float): Height in figure units.
    """
    import re

    body = re.sub(r"<\?xml[^>]*\?>", "", svg)
    view_box = re.search(r"viewBox=['\"]([^'\"]+)['\"]", body).group(1)
    body = re.sub(
        r"<svg[^>]*>",
        f'<svg x="{x}" y="{y}" width="{w}" height="{h}" viewBox="{view_box}" '
        f'preserveAspectRatio="xMidYMid meet">',
        body,
        count=1,
    )
    add(" ".join(body.split()))


def window(
    x: float, href: str, size: tuple[int, int] = (PIC_W, PIC_H), top: float = PIC_TOP
) -> None:
    """
    Place a transparent picture in a card's icon area, directly on the card's tint.

    Args:
        x (float): Card's left edge.
        href (str): Data URI of the picture.
        size (tuple[int, int]): Window width and height in figure px.
        top (float): Window's top, below the card's top edge.
    """
    width, height = size
    left = x + (W - width) / 2
    add(
        f'<image x="{left}" y="{Y + top}" width="{width}" height="{height}" href="{href}" '
        f'preserveAspectRatio="xMidYMid meet"/>'
    )


def card(x, y, w, h, title, bullets, accent):
    shape = f'x="{x}" y="{y}" width="{w}" height="{h}" rx="14"'
    clip = f"card-{x:.0f}-{y:.0f}"
    add(f'<clipPath id="{clip}"><rect {shape}/></clipPath>')
    add(f'<rect {shape} fill="{TINT[accent]}"/>')
    # The accent bar is clipped to the card, so it follows the rounded corners instead of
    # sticking out past them; the border is drawn last, over the bar's edge.
    add(f'<rect x="{x}" y="{y}" width="{w}" height="8" fill="{accent}" clip-path="url(#{clip})"/>')
    add(f'<rect {shape} fill="none" stroke="{accent}" stroke-width="2"/>')
    text(x + w / 2, y + 44, title, TITLE, INK, "bold", "middle")
    for i, line in enumerate(bullets):
        text(x + w / 2, y + 196 + 30 * i, line, BODY, INK_2, anchor="middle")


# ---------------------------------------------------------------- canvas
add(
    f'<svg xmlns="http://www.w3.org/2000/svg" width="180mm" height="{180 * HEIGHT / WIDTH:.1f}mm" '
    f'viewBox="0 0 {WIDTH} {HEIGHT}">'
)
add(
    f'<defs><marker id="head" viewBox="0 0 10 10" refX="8" refY="5" markerWidth="5" '
    f'markerHeight="5" orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" fill="{MUTED}"/>'
    f"</marker></defs>"
)
add(f'<rect width="{WIDTH}" height="{HEIGHT}" fill="{SURFACE}"/>')

# ---------------------------------------------------------------- panel a: pipeline
text(30, 40, "a", PANEL, INK, "bold")
text(70, 40, "Docking pipeline and the preparation choices it exposes", TITLE, INK, "bold")

W, H, GAP, Y = 222, 310, 43.6, 70
xs = [30 + i * (W + GAP) for i in range(6)]
cards = [
    ("Input", ["PoseBusters set", "428 complexes", "Crystal pose"], INK_2),
    ("Receptor", ["pH 7.4", "Cofactors removed", "No waters in set", "Held fixed"], INK_2),
    ("Ligand", ["Protonation", "(supplied or pH 7.4)", "Start conformer"], LIGAND),
    ("Search space", ["Box padding", "2 / 4 / 6 / 8 Å"], BOX),
    ("Engine", ["AutoDock Vina", "Smina", "Gnina", "rDock"], ENGINE),
    ("Evaluation", ["RMSD &lt; 2 Å", "PoseBusters", "validity"], INK_2),
]
if ARGS.v4:
    # Version 4 states the shown pose's RMSD, so the success criterion is concrete.
    rmsd = top_pose_rmsd(ARGS.complex)
    verdict = "✓" if rmsd < 2 else "✗"
    cards[5] = ("Evaluation", [f"RMSD {rmsd:.2f} Å {verdict}", "(success &lt; 2 Å)",
                               "PoseBusters validity"], INK_2)  # fmt: skip
    # The first line is left for the charge label, drawn in the ligand colour below.
    cards[2] = ("Ligand", ["", "Protonation", "(supplied or pH 7.4)", "Start conformer"], LIGAND)
    cards[3] = (
        "Search space",
        ["Box padding", "2 / 4 / 6 / 8 Å", "Baseline 4 Å", "Around crystal ligand"],
        BOX,
    )
for i, (x, (title, bullets, accent)) in enumerate(zip(xs, cards, strict=True), 1):
    card(x, Y, W, H, title, bullets, accent)
    if ARGS.v4:
        # Numbered badge on the card's top edge, clear of the title, so the flow reads at a
        # glance.
        add(
            f'<circle cx="{x + 26}" cy="{Y + 4}" r="15" fill="{accent}" stroke="#ffffff" '
            f'stroke-width="3"/>'
        )
        text(x + 26, Y + 11, str(i), 19, "#ffffff", "bold", "middle")
for x in xs[:-1]:
    arrow(x + W + 4, Y + 125, x + W + GAP - 6, Y + 125)

# Icons, centred in each card's upper half.
cy = Y + 120
c = [x + W / 2 for x in xs]
if ARGS.renders:
    # Version 3: a real complex. 3D renders for the receptor, search space, engine and
    # evaluation stages; the ligand's 2D structure for the ligand stage.
    for x, scene in zip(
        xs, ["input", "receptor", None, "box", "engine", "evaluation"], strict=True
    ):
        if scene:
            if ARGS.v4 and scene == "box":
                # Version 4: the box fills a larger window; the protein runs off its edges.
                size = (BOX_W, BOX_H)
                window(x, png_data(ARGS.renders / "box.png", size, BOX), size, BOX_TOP)
            elif ARGS.v4 and scene == "engine":
                # Version 4: the same window as the search space, so the poses fill it.
                size = (BOX_W, BOX_H)
                window(x, png_data(ARGS.renders / "engine.png", size), size, BOX_TOP)
            else:
                window(x, png_data(ARGS.renders / f"{scene}.png"))
    if ARGS.v4:
        # Version 4: the ligand card shows the change the protonation arm makes.
        drawing, charge = ligand_states_svg_data(ARGS.complex)
        inline_svg(drawing, xs[2] + (W - LIGAND_W) / 2, Y + PIC_TOP, LIGAND_W, PIC_H)
        # The charge is the first text line, so the card's lines match the receptor card's.
        text(xs[2] + W / 2, Y + 196, f"{charge:+d} at pH 7.4", BODY, LIGAND, "bold", "middle")
    else:
        inline_svg(
            ligand_svg_data(ARGS.complex), xs[2] + (W - PIC_W) / 2, Y + PIC_TOP, PIC_W, PIC_H
        )
else:
    protein(c[0], cy, MUTED, "#ececea")
    ligand(c[0] - 6, cy - 28, INK_2, "#ffffff", r=11)

    protein(c[1], cy, MUTED, "#ffffff")  # held fixed: drawn in neutral ink, not a stage colour
    for dx, dy in [(-44, -22), (-50, 18), (40, -26), (50, 14), (-8, 30), (24, 32)]:
        add(f'<circle cx="{c[1] + dx}" cy="{cy + dy}" r="5" fill="{MUTED}"/>')  # added H atoms
    text(c[1], cy + 6, "+H", 20, MUTED, "bold", "middle")

    ligand(c[2] - 14, cy, LIGAND, "#ffffff", r=22)
    add(
        f'<circle cx="{c[2] + 40}" cy="{cy - 22}" r="12" fill="#ffffff" '
        f'stroke="{LIGAND}" stroke-width="2.5"/>'
    )
    text(c[2] + 40, cy - 15, "±", 20, LIGAND, "bold", "middle")

    add(  # padded box (dashed) around the ligand's extent (solid)
        f'<rect x="{c[3] - 62}" y="{cy - 48}" width="124" height="96" fill="none" '
        f'stroke="{BOX}" stroke-width="2.5" stroke-dasharray="7 5"/>'
    )
    add(
        f'<rect x="{c[3] - 32}" y="{cy - 20}" width="64" height="40" fill="#ffffff" '
        f'stroke="{BOX}" stroke-width="2"/>'
    )
    ligand(c[3] - 14, cy, INK_2, r=10)
    arrow(c[3] + 34, cy, c[3] + 58, cy, BOX, 2)

    protein(c[4], cy + 14, ENGINE, "#ffffff")
    ligand(c[4] - 8, cy - 14, ENGINE, "#ffffff", r=11)
    add(  # sampling arrow: poses explored around the pocket
        f'<path d="M{c[4] - 58},{cy - 30} Q{c[4] - 36},{cy - 54} {c[4] - 14},{cy - 40}" '
        f'fill="none" stroke="{ENGINE}" stroke-width="2.5" marker-end="url(#head)"/>'
    )

    ligand(c[5] - 12, cy - 6, MUTED, r=20, dashed=True)  # crystal reference
    ligand(c[5] - 4, cy + 4, INK, r=20)  # docked pose

# Brackets naming what varies.
by = Y + H + 22
for (x0, x1), label, color in [
    ((xs[1], xs[1] + W), "Held fixed", MUTED),
    ((xs[2], xs[3] + W), "Varied one at a time", INK_2),
    ((xs[4], xs[4] + W), "Engine", ENGINE),
]:
    add(
        f'<path d="M{x0 + 6},{by} v12 H{x1 - 6} v-12" fill="none" '
        f'stroke="{color}" stroke-width="2"/>'
    )
    text((x0 + x1) / 2, by + 42, label, BODY, color, "bold", "middle")

# ---------------------------------------------------------------- panel b: design
PY = 520
text(30, PY, "b", PANEL, INK, "bold")
text(70, PY, "One factor at a time, paired by complex", TITLE, INK, "bold")

arms = [
    ("Baseline (4 Å box)", INK_2),
    ("Box padding 2 Å", BOX),
    ("Box padding 6 Å", BOX),
    ("Box padding 8 Å", BOX),
    ("Ligand protonation", LIGAND),
    ("Start conformer", LIGAND),
]
engines = ["Vina", "Smina", "Gnina", "rDock"]
MX, MY, CW, CH, PITCH_X, PITCH_Y = 290, 590, 58, 34, 70, 46
for j, name in enumerate(engines):
    text(MX + j * PITCH_X + CW / 2, MY - 14, name, 20, ENGINE, "bold", "middle")  # 6.4 pt; fits
for i, (arm, color) in enumerate(arms):
    y = MY + i * PITCH_Y
    add(f'<rect x="40" y="{y + 9}" width="8" height="16" rx="2" fill="{color}"/>')
    text(58, y + 25, arm, BODY, INK if i == 0 else INK_2, "bold" if i == 0 else "normal")
    for j in range(len(engines)):
        fill = "#c3c2b7" if i == 0 else GRID
        add(
            f'<rect x="{MX + j * PITCH_X}" y="{y}" width="{CW}" height="{CH}" '
            f'rx="5" fill="{fill}"/>'
        )
# Highlight one comparison: ligand protonation vs baseline under Vina.
for i in (0, 4):
    add(
        f'<rect x="{MX - 3}" y="{MY + i * PITCH_Y - 3}" width="{CW + 6}" height="{CH + 6}" rx="7" '
        f'fill="none" stroke="{LIGAND}" stroke-width="3"/>'
    )
# Below both the arm list and the 8 paired-outcome rows drawn next to it (36 px apart).
legend_y = max(MY + len(arms) * PITCH_Y, MY + 30 + 7 * 36 + 10) + 30
add(f'<rect x="40" y="{legend_y - 20}" width="30" height="20" rx="4" fill="{GRID}"/>')
text(80, legend_y - 3, "one run: 428 complexes", BODY, INK_2)

# Paired outcomes for the highlighted comparison.
PX0, PX1, DY0 = 650, 740, MY + 30
text(PX0, MY - 14, "Base", 22, INK_2, "bold", "middle")
text(PX1, MY - 14, "Arm", 22, LIGAND, "bold", "middle")
outcomes = [(1, 1), (1, 0), (0, 0), (1, 1), (0, 1), (1, 1), (0, 0), (0, 1)]
for k, (a, b) in enumerate(outcomes):
    y = DY0 + k * 36
    discordant = a != b
    add(
        f'<line x1="{PX0 + 11}" y1="{y}" x2="{PX1 - 11}" y2="{y}" '
        f'stroke="{LIGAND if discordant else GRID}" stroke-width="{3 if discordant else 2}"/>'
    )
    for x, hit in ((PX0, a), (PX1, b)):
        fill = INK_2 if hit else "#ffffff"
        add(f'<circle cx="{x}" cy="{y}" r="10" fill="{fill}" stroke="{INK_2}" stroke-width="2"/>')
# Legend: symbols centred on one column, labels starting on another.
KEY_X, LABEL_X, ROW = 625, 650, 34
key_rows = [
    (f'<circle cx="{KEY_X}" cy="{{y}}" r="10" fill="{INK_2}" stroke="{INK_2}" stroke-width="2"/>',
     "success"),
    (f'<circle cx="{KEY_X}" cy="{{y}}" r="10" fill="#ffffff" stroke="{INK_2}" stroke-width="2"/>',
     "miss"),
    (f'<line x1="{KEY_X - 15}" y1="{{y}}" x2="{KEY_X + 15}" y2="{{y}}" stroke="{LIGAND}" '
     f'stroke-width="3"/>', "outcome changed"),
]  # fmt: skip
for k, (symbol, label) in enumerate(key_rows):
    y = legend_y - 10 + k * ROW
    add(symbol.format(y=y))
    text(LABEL_X, y + 7, label, BODY, INK_2)

# ---------------------------------------------------------------- panel c: question
CX0 = 880
text(CX0 - 40, PY, "c", PANEL, INK, "bold")
text(CX0, PY, "Which matters more: preparation or engine?", TITLE, INK, "bold")

AX0, AX1 = CX0 + 20, 1550


def dot_row(y, label, dots, colors, delta_label, delta_color):
    text(AX0, y - 46, label, BODY, INK_2)
    add(f'<line x1="{AX0}" y1="{y}" x2="{AX1}" y2="{y}" stroke="{GRID}" stroke-width="2"/>')
    for x, color in zip(dots, colors, strict=True):
        add(f'<circle cx="{x}" cy="{y}" r="11" fill="{color}" stroke="#ffffff" stroke-width="2"/>')
    lo, hi = min(dots), max(dots)
    add(
        f'<path d="M{lo},{y + 22} v10 H{hi} v-10" fill="none" '
        f'stroke="{delta_color}" stroke-width="2"/>'
    )
    text((lo + hi) / 2, y + 60, delta_label, BODY, delta_color, "bold", "middle")


dot_row(
    650,
    "Engines, at a fixed preparation",
    [1080, 1160, 1215, 1290],
    [ENGINE] * 4,
    "Δ engine",
    ENGINE,
)
dot_row(
    800,
    "Preparation choices, at a fixed engine",
    [1040, 1120, 1190, 1260, 1330],
    [BOX, BOX, LIGAND, LIGAND, BOX],
    "Δ preparation",
    INK_2,
)
arrow(AX0, 905, AX1, 905, MUTED)
text((AX0 + AX1) / 2, 940, "Docking success rate", BODY, INK_2, anchor="middle")

add("</svg>")

default = "concept_v4.svg" if ARGS.v4 else "concept_v3.svg" if ARGS.renders else "concept_v2.svg"
target = ARGS.out or Path("docs/figures") / default
target.parent.mkdir(parents=True, exist_ok=True)
target.write_text("\n".join(out) + "\n", encoding="utf-8", newline="\n")  # same on every OS
print(target)
