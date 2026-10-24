"""Draw Figure 1 version 5 (study design) as an illustrative vector SVG.

Usage: python scripts/make_concept_figure_v5.py [OUT.svg]   (default docs/figures/concept_v5.svg)

Version 5 replaces the renders and line drawings of versions 1-4 with one illustrative style:
six large panels in a 3 x 2 grid, a shaded protein with a pocket, stylised ligands and a 3D
search box, one short caption per panel, and fixed colour roles (protein slate blue, ligand
orange, search box green, docked poses and engines purple, crystal pose amber).

The figure is 180 mm wide: 1 px = 0.32 pt, so 22 px text prints at 7 pt.
"""

import math
import sys
from pathlib import Path

WIDTH, HEIGHT = 1600, 1040
FONT = "Arial, Helvetica, sans-serif"

INK, INK_2, MUTED, GRID = "#1d1d1b", "#52514e", "#898781", "#e1e0d9"
PANEL_BG, PANEL_LINE = "#f7f8fa", "#d9dde3"
PROTEIN, PROTEIN_DARK = "#5b84b1", "#33587f"
LIGAND, BOX, ENGINE, CRYSTAL = "#eb6834", "#1baf7a", "#6250d6", "#f0a830"
SUCCESS = "#2e9e4f"

out: list[str] = []


def add(element: str) -> None:
    """
    Append one SVG element.

    Args:
        element (str): SVG markup.
    """
    out.append(element)


def text(x, y, s, size=22, color=INK_2, weight="normal", anchor="middle", style="normal"):
    """
    Add a line of text.

    Args:
        x (float): Anchor x.
        y (float): Baseline y.
        s (str): Text (XML-escaped by the caller).
        size (int): Font size in px.
        color (str): Fill colour.
        weight (str): Font weight.
        anchor (str): text-anchor.
        style (str): font-style.
    """
    add(
        f'<text x="{x:.1f}" y="{y:.1f}" font-family="{FONT}" font-size="{size}" '
        f'font-weight="{weight}" font-style="{style}" fill="{color}" '
        f'text-anchor="{anchor}">{s}</text>'
    )


def panel(x, y, w, h, number, title, accent, tag=None):
    """
    Draw a panel: rounded card, numbered badge and title, optional tag chip at the right.

    Args:
        x (float): Left edge.
        y (float): Top edge.
        w (float): Width.
        h (float): Height.
        number (int): Panel number shown in the badge.
        title (str): Panel title.
        accent (str): Badge colour.
        tag (str | None): Small grey chip after the title, e.g. "held fixed".
    """
    add(
        f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="22" fill="{PANEL_BG}" '
        f'stroke="{PANEL_LINE}" stroke-width="2"/>'
    )
    add(f'<circle cx="{x + 42}" cy="{y + 44}" r="20" fill="{accent}"/>')
    text(x + 42, y + 52, str(number), 24, "#ffffff", "bold")
    text(x + 74, y + 54, title, 30, INK, "bold", "start")
    if tag:
        chip(x + w - 28, y + 44, tag, MUTED, "#ffffff", anchor="end")


def chip(x, y, label, color, fill, anchor="middle", size=19):
    """
    Draw a rounded label chip centred vertically on y.

    Args:
        x (float): Anchor x (centre, or right edge with anchor="end").
        y (float): Vertical centre.
        label (str): Chip text.
        color (str): Text and outline colour.
        fill (str): Background.
        anchor (str): "middle" or "end".
        size (int): Font size.

    Returns:
        float: Chip width.
    """
    w = 0.56 * size * len(label) + 24
    left = x - w / 2 if anchor == "middle" else x - w
    add(
        f'<rect x="{left:.1f}" y="{y - 16}" width="{w:.1f}" height="32" rx="16" fill="{fill}" '
        f'stroke="{color}" stroke-width="1.6"/>'
    )
    text(left + w / 2, y + 7, label, size, color, "bold")
    return w


def bean(cx, cy, s, pocket=True):
    """
    Draw the protein: a shaded rounded body with a binding pocket in its top edge.

    Args:
        cx (float): Centre x.
        cy (float): Centre y.
        s (float): Scale (1 = 124 x 80 px).
        pocket (bool): Notch a pocket into the top edge.
    """
    notch = "L-16,-14 Q0,-2 16,-14 L20,-38" if pocket else "C-6,-40 6,-40 20,-38"
    d = (
        f"M-62,8 C-62,-30 -34,-40 -20,-38 {notch} C42,-40 62,-24 62,6 "
        "C62,36 30,42 0,40 C-36,42 -62,34 -62,8 Z"
    )
    add(
        f'<g transform="translate({cx},{cy}) scale({s})">'
        f'<path d="{d}" fill="url(#protein)" stroke="{PROTEIN_DARK}" stroke-width="{2.4 / s:.2f}"/>'
        f'<ellipse cx="-34" cy="-14" rx="14" ry="7" fill="#ffffff" opacity="0.35" '
        f'transform="rotate(-25 -34 -14)"/></g>'
    )


def hexagon(cx, cy, r, rotation=0):
    """
    Return the points of a regular hexagon.

    Args:
        cx (float): Centre x.
        cy (float): Centre y.
        r (float): Circumradius.
        rotation (float): Rotation in degrees.

    Returns:
        str: Points attribute.
    """
    return " ".join(
        f"{cx + r * math.cos(math.radians(a + rotation)):.1f},"
        f"{cy + r * math.sin(math.radians(a + rotation)):.1f}"
        for a in range(30, 390, 60)
    )


def molecule(cx, cy, s, color, width=4.5, dashed=False, fill="#ffffff", angle=0, charge=False):
    """
    Draw the stylised ligand: two rings joined by a short linker, with a side group.

    Args:
        cx (float): Centre x.
        cy (float): Centre y.
        s (float): Scale (1 = about 110 x 50 px).
        color (str): Stroke colour.
        width (float): Stroke width at scale 1.
        dashed (bool): Dashed outline (e.g. an alternative conformer).
        fill (str): Ring fill.
        angle (float): Rotation in degrees.
        charge (bool): Add a "+" charge badge on the side group.
    """
    dash = ' stroke-dasharray="7 5"' if dashed else ""
    stroke = (
        f'stroke="{color}" stroke-width="{width}" stroke-linejoin="round" '
        f'stroke-linecap="round"{dash}'
    )
    parts = [
        f'<polygon points="{hexagon(-34, 4, 20)}" fill="{fill}" {stroke}/>',
        f'<polyline points="-15,-2 2,-12 20,-4" fill="none" {stroke}/>',
        f'<polygon points="{hexagon(38, 2, 20)}" fill="{fill}" {stroke}/>',
        f'<polyline points="-51,15 -60,32 -76,34" fill="none" {stroke}/>',
    ]
    if charge:
        parts.append(
            f'<circle cx="-82" cy="36" r="12" fill="#ffffff" stroke="{color}" stroke-width="2.5"/>'
        )
        parts.append(
            f'<text x="-82" y="43" font-family="{FONT}" font-size="22" font-weight="bold" '
            f'fill="{color}" text-anchor="middle">+</text>'
        )
    add(f'<g transform="translate({cx},{cy}) rotate({angle}) scale({s})">{"".join(parts)}</g>')


def cube(cx, cy, size, color, width=3, dash=None, grid=0, depth=0.42):
    """
    Draw a box in oblique projection, optionally with a grid on its front face.

    Args:
        cx (float): Centre of the front face, x.
        cy (float): Centre of the front face, y.
        size (float): Edge length of the front face.
        color (str): Line colour.
        width (float): Line width.
        dash (str | None): stroke-dasharray.
        grid (int): Grid divisions on the front face (0 for none).
        depth (float): Receding edge length as a fraction of size.
    """
    h = size / 2
    dx, dy = size * depth * 0.8, -size * depth * 0.6
    style = f'fill="none" stroke="{color}" stroke-width="{width}" stroke-linejoin="round"'
    if dash:
        style += f' stroke-dasharray="{dash}"'
    x0, y0, x1, y1 = cx - h, cy - h, cx + h, cy + h
    back = f"M{x0 + dx},{y0 + dy} H{x1 + dx} V{y1 + dy}"
    add(f'<path d="{back}" {style} opacity="0.55"/>')
    for x, y in ((x0, y0), (x1, y0), (x1, y1)):
        add(f'<line x1="{x}" y1="{y}" x2="{x + dx}" y2="{y + dy}" {style}/>')
    add(f'<rect x="{x0}" y="{y0}" width="{size}" height="{size}" {style}/>')
    for i in range(1, grid):
        t = x0 + size * i / grid
        u = y0 + size * i / grid
        add(
            f'<line x1="{t}" y1="{y0}" x2="{t}" y2="{y1}" stroke="{color}" '
            f'stroke-width="1" opacity="0.45"/>'
        )
        add(
            f'<line x1="{x0}" y1="{u}" x2="{x1}" y2="{u}" stroke="{color}" '
            f'stroke-width="1" opacity="0.45"/>'
        )


def arrow(x1, y1, x2, y2, color=MUTED, width=4):
    """
    Draw an arrow with the shared arrowhead marker.

    Args:
        x1 (float): Start x.
        y1 (float): Start y.
        x2 (float): End x.
        y2 (float): End y.
        color (str): Colour (the marker is grey).
        width (float): Line width.
    """
    add(
        f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="{color}" stroke-width="{width}" '
        f'stroke-linecap="round" marker-end="url(#head)"/>'
    )


def curved_arrow(x1, y1, x2, y2, bend, color=ENGINE, width=3):
    """
    Draw a curved arrow (sampling around the pocket).

    Args:
        x1 (float): Start x.
        y1 (float): Start y.
        x2 (float): End x.
        y2 (float): End y.
        bend (float): Control-point offset perpendicular to the chord.
        color (str): Colour.
        width (float): Line width.
    """
    mx, my = (x1 + x2) / 2, (y1 + y2) / 2
    length = math.hypot(x2 - x1, y2 - y1) or 1
    nx, ny = -(y2 - y1) / length, (x2 - x1) / length
    add(
        f'<path d="M{x1},{y1} Q{mx + nx * bend:.1f},{my + ny * bend:.1f} {x2},{y2}" fill="none" '
        f'stroke="{color}" stroke-width="{width}" stroke-linecap="round" marker-end="url(#headp)"/>'
    )


# ---------------------------------------------------------------- canvas
add(
    f'<svg xmlns="http://www.w3.org/2000/svg" width="180mm" height="{180 * HEIGHT / WIDTH:.1f}mm" '
    f'viewBox="0 0 {WIDTH} {HEIGHT}">'
)
add(
    "<defs>"
    '<radialGradient id="protein" cx="0.35" cy="0.3" r="0.85">'
    '<stop offset="0" stop-color="#a9c5e6"/><stop offset="0.55" stop-color="#6f96c3"/>'
    f'<stop offset="1" stop-color="{PROTEIN}"/></radialGradient>'
    '<marker id="head" viewBox="0 0 10 10" refX="7" refY="5" markerWidth="4" markerHeight="4" '
    f'orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" fill="{MUTED}"/></marker>'
    '<marker id="headp" viewBox="0 0 10 10" refX="7" refY="5" markerWidth="4.5" markerHeight="4.5" '
    f'orient="auto"><path d="M0,0 L10,5 L0,10 z" fill="{ENGINE}"/></marker>'
    "</defs>"
)
add(f'<rect width="{WIDTH}" height="{HEIGHT}" fill="#ffffff"/>')

M, GAP, PW, PH = 16, 56, 484, 476
xs = [M + i * (PW + GAP) for i in range(3)]
Y1, Y2 = 16, 16 + PH + 56

# ---------------------------------------------------------------- 1 receptor (held fixed)
x, y = xs[0], Y1
panel(x, y, PW, PH, 1, "Receptor", MUTED, tag="held fixed")
cx, cy = x + PW / 2, y + 225
bean(cx, cy, 2.6)
for ax, ay in [(-150, -40), (-160, 30), (-90, 98), (20, 104), (120, 84), (158, 10), (140, -70)]:
    add(
        f'<circle cx="{cx + ax}" cy="{cy + ay}" r="9" fill="#ffffff" '
        f'stroke="{PROTEIN_DARK}" stroke-width="2.4"/>'
    )
    text(cx + ax, cy + ay + 5, "H", 13, PROTEIN_DARK, "bold")
text(cx, y + 410, "Hydrogens added at pH 7.4,", 22)
text(cx, y + 440, "cofactors removed (no waters in the set)", 22)

# ---------------------------------------------------------------- 2 ligand
x, y = xs[1], Y1
panel(x, y, PW, PH, 2, "Ligand", LIGAND)
# Protonation: supplied state -> dominant state at pH 7.4 (charge badge on the side group).
molecule(x + 140, y + 150, 1.1, LIGAND)
arrow(x + 228, y + 150, x + 266, y + 150)
molecule(x + 370, y + 150, 1.1, LIGAND, charge=True)
text(x + 140, y + 232, "supplied", 20, INK_2)
text(x + 362, y + 232, "pH 7.4", 20, INK_2)
add(
    f'<line x1="{x + 40}" y1="{y + 262}" x2="{x + PW - 40}" y2="{y + 262}" '
    f'stroke="{GRID}" stroke-width="2"/>'
)
# Start conformer: supplied, or regenerated with RDKit.
molecule(x + 140, y + 318, 1.0, LIGAND)
text(x + PW / 2 + 2, y + 326, "or", 22, MUTED, "bold")
molecule(x + 352, y + 318, 1.0, LIGAND, dashed=True, angle=-18)
text(x + 140, y + 384, "supplied", 20, INK_2)
text(x + 352, y + 384, "RDKit", 20, INK_2)
text(x + PW / 2, y + 440, "Varied: protonation, start conformer", 22)

# ---------------------------------------------------------------- 3 search space
x, y = xs[2], Y1
panel(x, y, PW, PH, 3, "Search space", BOX)
cx, cy = x + PW / 2 - 20, y + 215
bean(cx + 10, cy + 62, 2.0)
bx, by = cx, cy  # front-face centre of the boxes, over the pocket
# Nested boxes, one per padding; the baseline (4 Å) solid with a grid. A key below names
# them, so the drawing carries no labels.
for pad, size in [(8, 186), (6, 154), (4, 122), (2, 90)]:
    if pad == 4:
        cube(bx, by, size, BOX, width=3.5, grid=5, depth=0.32)
    else:
        cube(bx, by, size, BOX, width=2, dash="6 6", depth=0.32)
molecule(bx, by + 6, 0.6, LIGAND, width=6)
ky = y + 368
add(f'<line x1="{x + 92}" y1="{ky}" x2="{x + 126}" y2="{ky}" stroke="{BOX}" stroke-width="3.5"/>')
text(x + 134, ky + 7, "4 Å (baseline)", 19, INK_2, anchor="start")
add(
    f'<line x1="{x + 290}" y1="{ky}" x2="{x + 324}" y2="{ky}" stroke="{BOX}" '
    'stroke-width="2" stroke-dasharray="6 6"/>'
)
text(x + 332, ky + 7, "2, 6, 8 Å", 19, INK_2, anchor="start")
text(x + PW / 2, y + 430, "Box around the crystal ligand", 22)

# ---------------------------------------------------------------- 4 docking
x, y = xs[0], Y2
panel(x, y, PW, PH, 4, "Docking", ENGINE)
cx, cy = x + PW / 2, y + 250
bean(cx, cy + 20, 2.4)
for dx, dy, a in [(-120, -120, -30), (120, -125, 25), (-150, -40, 60)]:
    molecule(cx + dx, cy + dy, 0.62, "#a99fe8", width=5, angle=a)
molecule(cx, cy - 62, 0.7, ENGINE, width=6)
curved_arrow(cx - 95, cy - 150, cx - 30, cy - 92, -30)
curved_arrow(cx + 95, cy - 152, cx + 32, cy - 92, 30)
curved_arrow(cx - 152, cy + 4, cx - 72, cy - 48, 30)
chips = ["Vina", "Smina", "Gnina", "rDock"]
widths = [0.56 * 19 * len(c) + 24 for c in chips]
left = cx - (sum(widths) + 12 * 3) / 2
for c, w in zip(chips, widths, strict=True):
    chip(left + w / 2, y + 392, c, ENGINE, "#f1effc")
    left += w + 12
text(cx, y + 444, "Four engines, nine poses each", 22)

# ---------------------------------------------------------------- 5 evaluation
x, y = xs[1], Y2
panel(x, y, PW, PH, 5, "Evaluation", INK_2)
cx, cy = x + PW / 2, y + 190
molecule(cx, cy, 1.7, CRYSTAL, width=11, fill="#fdf1dc")
molecule(cx + 5, cy + 3, 1.7, ENGINE, width=4, fill="none", angle=3)
text(cx, y + 290, "RMSD &lt; 2 Å", 26, INK, "bold")
add(
    f'<path d="M{cx + 92},{y + 280} l9,10 l18,-22" fill="none" stroke="{SUCCESS}" '
    'stroke-width="5" stroke-linecap="round" stroke-linejoin="round"/>'
)
checks = ["bonds", "angles", "clashes", "planarity"]
widths = [0.56 * 19 * len(c) + 24 for c in checks]
left = cx - (sum(widths) + 10 * 3) / 2
for c, w in zip(checks, widths, strict=True):
    chip(left + w / 2, y + 352, c, SUCCESS, "#eef8f1")
    left += w + 10
text(cx, y + 410, "Top pose vs crystal pose (amber);", 22)
text(cx, y + 440, "PB-valid if it passes PoseBusters", 22)

# ---------------------------------------------------------------- 6 paired comparison
x, y = xs[2], Y2
panel(x, y, PW, PH, 6, "Compare, by complex", LIGAND)
bx, ax_, top = x + 70, x + 160, y + 120
text(bx, top - 14, "base", 19, INK_2, "bold")
text(ax_, top - 14, "arm", 19, LIGAND, "bold")
for k, (a, b) in enumerate([(1, 1), (1, 0), (0, 0), (0, 1), (1, 1), (0, 1)]):
    yy = top + 12 + k * 40
    changed = a != b
    add(
        f'<line x1="{bx + 12}" y1="{yy}" x2="{ax_ - 12}" y2="{yy}" '
        f'stroke="{LIGAND if changed else GRID}" stroke-width="{3.5 if changed else 2}"/>'
    )
    for px, hit in ((bx, a), (ax_, b)):
        fill = INK_2 if hit else "#ffffff"
        add(
            f'<circle cx="{px}" cy="{yy}" r="11" fill="{fill}" stroke="{INK_2}" '
            f'stroke-width="2.2"/>'
        )
# Δ engine vs Δ preparation, schematic.
cx0, cx1 = x + 225, x + PW - 30
for row, (label, dots, colors, bracket_color) in enumerate(
    [
        ("Δ engine", [cx0 + 30, cx0 + 105, cx0 + 140, cx0 + 215], [ENGINE] * 4, ENGINE),
        ("Δ preparation", [cx0 + 75, cx0 + 100, cx0 + 125, cx0 + 150],
         [BOX, LIGAND, BOX, LIGAND], INK_2),
    ]
):  # fmt: skip
    yy = y + 160 + row * 120
    add(f'<line x1="{cx0}" y1="{yy}" x2="{cx1}" y2="{yy}" stroke="{GRID}" stroke-width="2"/>')
    for px, col in zip(dots, colors, strict=True):
        add(f'<circle cx="{px}" cy="{yy}" r="10" fill="{col}" stroke="#ffffff" stroke-width="2"/>')
    lo, hi = min(dots), max(dots)
    add(
        f'<path d="M{lo},{yy + 20} v8 H{hi} v-8" fill="none" stroke="{bracket_color}" '
        f'stroke-width="2.4"/>'
    )
    text((lo + hi) / 2, yy + 54, label, 20, bracket_color, "bold")
# Legend for the paired dots.
ly = y + 382
add(f'<circle cx="{x + 44}" cy="{ly}" r="9" fill="{INK_2}"/>')
text(x + 60, ly + 7, "success", 19, INK_2, anchor="start")
add(f'<circle cx="{x + 160}" cy="{ly}" r="9" fill="#ffffff" stroke="{INK_2}" stroke-width="2"/>')
text(x + 176, ly + 7, "miss", 19, INK_2, anchor="start")
add(
    f'<line x1="{x + 250}" y1="{ly}" x2="{x + 280}" y2="{ly}" stroke="{LIGAND}" '
    f'stroke-width="3.5"/>'
)
text(x + 290, ly + 7, "outcome changed", 19, INK_2, anchor="start")
text(x + PW / 2, y + 440, "Each arm vs its baseline, paired", 22)

# ---------------------------------------------------------------- flow arrows
mid1 = Y1 + PH / 2
for i in range(2):
    arrow(xs[i] + PW + 10, mid1, xs[i + 1] - 12, mid1)
mid2 = Y2 + PH / 2
for i in range(2):
    arrow(xs[i] + PW + 10, mid2, xs[i + 1] - 12, mid2)
# From panel 3 down to panel 4: under the top row, right to left.
ty = Y1 + PH + 28
add(
    f'<path d="M{xs[2] + PW / 2},{Y1 + PH + 6} V{ty} H{xs[0] + PW / 2} V{Y2 - 10}" fill="none" '
    f'stroke="{MUTED}" stroke-width="4" stroke-linejoin="round" marker-end="url(#head)"/>'
)

add("</svg>")
target = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("docs/figures/concept_v5.svg")
target.parent.mkdir(parents=True, exist_ok=True)
target.write_text("\n".join(out) + "\n", encoding="utf-8", newline="\n")  # same on every OS
print(target)
