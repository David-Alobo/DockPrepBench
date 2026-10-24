"""Draw the study-design figure (Figure 1) as SVG.

Usage: python scripts/make_concept_figure.py [OUT.svg]

Version 2 (default output docs/figures/concept_v2.svg) shows the study as run: the receptor
is held fixed, and the arms vary the search box and the ligand. Version 1 (concept.svg, the
original plan with receptor, cofactor and water arms) is kept for reference.

The figure is 180 mm wide (double column). At that size 1 px = 0.32 pt, so body text
(22 px) prints at 7 pt. Convert to PDF/PNG with e.g. `rsvg-convert`.
"""

import math
import sys
from pathlib import Path

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


def card(x, y, w, h, title, bullets, accent):
    add(
        f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="14" fill="{TINT[accent]}" '
        f'stroke="{accent}" stroke-width="2"/>'
    )
    add(f'<rect x="{x}" y="{y}" width="{w}" height="8" rx="4" fill="{accent}"/>')
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
for x, (title, bullets, accent) in zip(xs, cards, strict=True):
    card(x, Y, W, H, title, bullets, accent)
for x in xs[:-1]:
    arrow(x + W + 4, Y + 125, x + W + GAP - 6, Y + 125)

# Icons, centred in each card's upper half.
cy = Y + 120
c = [x + W / 2 for x in xs]
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
    f'<path d="M{c[4] - 58},{cy - 30} Q{c[4] - 36},{cy - 54} {c[4] - 14},{cy - 40}" fill="none" '
    f'stroke="{ENGINE}" stroke-width="2.5" marker-end="url(#head)"/>'
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
    text(MX + j * PITCH_X + CW / 2, MY - 14, name, 20, ENGINE, "bold", "middle")
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
text(PX0, MY - 14, "Base", 20, INK_2, "bold", "middle")
text(PX1, MY - 14, "Arm", 20, LIGAND, "bold", "middle")
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
add(f'<circle cx="610" cy="{legend_y - 10}" r="10" fill="{INK_2}"/>')
text(628, legend_y - 3, "success", BODY, INK_2)
add(
    f'<line x1="610" y1="{legend_y + 26}" x2="640" y2="{legend_y + 26}" '
    f'stroke="{LIGAND}" stroke-width="3"/>'
)
text(650, legend_y + 33, "changed", BODY, INK_2)

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
text(
    (AX0 + AX1) / 2,
    985,
    "Schematic: positions are illustrative, not results.",
    18,
    MUTED,
    anchor="middle",
    style="italic",
)

add("</svg>")

target = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("docs/figures/concept_v2.svg")
target.parent.mkdir(parents=True, exist_ok=True)
target.write_text("\n".join(out) + "\n", encoding="utf-8")
print(target)
