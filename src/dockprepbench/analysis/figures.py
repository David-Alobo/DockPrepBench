"""Publication figures comparing benchmark runs."""

from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib.figure import Figure

from dockprepbench.analysis.stats import success_rate, successes
from dockprepbench.evaluation.rmsd import SUCCESS_THRESHOLD

# Categorical slots in fixed order (validated light-mode palette): a run keeps its slot,
# and so its colour, in every figure. Never cycled: more runs than slots is an error.
SERIES = ("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948")
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_SECONDARY = "#52514e"
INK_MUTED = "#898781"
GRID = "#e1e0d9"
BASELINE = "#c3c2b7"
MAX_DIRECT_LABELS = 4  # above this, rely on the legend alone to avoid collisions
FORMATS = ("png", "pdf")


def summary_table(runs: dict[str, pd.DataFrame]) -> pd.DataFrame:
    """
    Tabulate top-1 and PB-valid success with 95% Wilson intervals for each run.

    Args:
        runs (dict[str, pd.DataFrame]): Run name to `load_run` output, in display order.

    Returns:
        pd.DataFrame: One row per run, indexed by name, rates as fractions.
    """
    rows = []
    for name, run in runs.items():
        top1 = success_rate(successes(run))
        valid = success_rate(successes(run, pb_valid=True))
        rows.append(
            {
                "run": name,
                "n": top1.total,
                "failed": int((run["error"] != "").sum()),
                "top1": top1.value,
                "top1_low": top1.low,
                "top1_high": top1.high,
                "pb_valid": valid.value,
                "pb_valid_low": valid.low,
                "pb_valid_high": valid.high,
            }
        )
    return pd.DataFrame(rows).set_index("run")


def plot_success_rates(table: pd.DataFrame, path_stem: Path) -> list[Path]:
    """
    Draw a forest plot of top-1 and PB-valid success per run, with 95% intervals.

    The first run is the reference: a dashed line marks its top-1 rate.

    Args:
        table (pd.DataFrame): Output of `summary_table`.
        path_stem (Path): Output path without extension; one file per format is written.

    Returns:
        list[Path]: Files written.
    """
    fig = Figure(figsize=(6.4, 0.55 * len(table) + 1.4), facecolor=SURFACE)
    ax = fig.add_subplot()
    _style_axes(ax)

    rows = np.arange(len(table))[::-1]  # first run at the top
    metrics = [("top1", "Top-1 (RMSD < 2 Å)", SERIES[0]), ("pb_valid", "PB-valid", SERIES[1])]
    for offset, (key, label, colour) in zip((0.14, -0.14), metrics, strict=True):
        y = rows + offset
        ax.hlines(y, 100 * table[f"{key}_low"], 100 * table[f"{key}_high"], color=colour, lw=2)
        ax.scatter(100 * table[key], y, s=64, color=colour, edgecolor=SURFACE, lw=2, zorder=3)
        ax.plot([], [], marker="o", ms=8, lw=2, color=colour, label=label)  # legend handle
        # Direct label on the top row only, so identity never rests on colour alone.
        _label(ax, label, 100 * table[f"{key}_high"].iloc[0], y[0])

    ax.axvline(100 * table["top1"].iloc[0], color=INK_MUTED, lw=1, ls=(0, (4, 3)), zorder=1)
    ax.set_yticks(rows, table.index)
    ax.set_ylim(-0.6, len(table) - 0.4)
    ax.set_xlim(0, 100)
    ax.set_xlabel("Complexes docked successfully (%)", color=INK_SECONDARY)
    ax.grid(axis="x", color=GRID, lw=0.8)
    _legend_above(ax, ncol=2)
    return _save(fig, path_stem)


def plot_rmsd_cdf(
    runs: dict[str, pd.DataFrame], path_stem: Path, max_rmsd: float = 10.0
) -> list[Path]:
    """
    Plot the cumulative fraction of complexes whose top pose is within each RMSD.

    Failed complexes never reach the curve, so every curve plateaus at its success
    ceiling rather than at 100%.

    Args:
        runs (dict[str, pd.DataFrame]): Run name to `load_run` output, in display order.
        path_stem (Path): Output path without extension; one file per format is written.
        max_rmsd (float): Upper end of the RMSD axis in Å.

    Returns:
        list[Path]: Files written.

    Raises:
        ValueError: If there are more runs than categorical colour slots.
    """
    if len(runs) > len(SERIES):
        raise ValueError(f"At most {len(SERIES)} runs per figure; facet or fold the rest")
    fig = Figure(figsize=(6.4, 4.2), facecolor=SURFACE)
    ax = fig.add_subplot()
    _style_axes(ax)

    grid = np.linspace(0, max_rmsd, 501)
    ends = {}
    for (name, run), colour in zip(runs.items(), SERIES, strict=False):
        rmsds = run["top1_rmsd"].to_numpy(dtype=float)
        fraction = [100 * np.mean(rmsds <= x) for x in grid]  # NaN compares False: a miss
        ax.plot(grid, fraction, color=colour, lw=2, label=name, drawstyle="steps-post")
        ends[name] = fraction[-1]

    if len(runs) <= MAX_DIRECT_LABELS:
        # Curves converge at the right edge, so spread the end labels apart.
        for name, (y_end, y_text) in zip(ends, _spread(list(ends.values()), 5.0), strict=True):
            ax.annotate(
                name,
                (max_rmsd, y_end),
                xytext=(max_rmsd + 0.3, y_text),
                va="center",
                fontsize=8,
                color=INK_SECONDARY,
                annotation_clip=False,
                arrowprops={"arrowstyle": "-", "color": GRID, "lw": 0.8},
            )

    ax.axvline(SUCCESS_THRESHOLD, color=INK_MUTED, lw=1, ls=(0, (4, 3)))
    _label(ax, f"{SUCCESS_THRESHOLD:g} Å", SUCCESS_THRESHOLD, 3, color=INK_MUTED)
    ax.set_xlim(0, max_rmsd)
    ax.set_ylim(0, 100)
    ax.set_xlabel("Top-1 pose RMSD (Å)", color=INK_SECONDARY)
    ax.set_ylabel("Complexes within RMSD (%)", color=INK_SECONDARY)
    ax.grid(color=GRID, lw=0.8)
    _legend_above(ax, ncol=min(len(runs), 4))
    return _save(fig, path_stem)


def _spread(values: list[float], min_gap: float) -> list[tuple[float, float]]:
    """
    Nudge label positions apart so neighbours are at least `min_gap` apart.

    Args:
        values (list[float]): Anchor positions, in data units.
        min_gap (float): Minimum distance between label positions.

    Returns:
        list[tuple[float, float]]: (anchor, label position) in the input order.
    """
    order = sorted(range(len(values)), key=values.__getitem__)
    placed = [values[i] for i in order]
    for k in range(1, len(placed)):
        placed[k] = max(placed[k], placed[k - 1] + min_gap)
    # Re-centre the group so it does not drift upwards as a whole.
    shift = (sum(placed) - sum(values[i] for i in order)) / len(placed)
    result = [0.0] * len(values)
    for k, i in enumerate(order):
        result[i] = placed[k] - shift
    return list(zip(values, result, strict=True))


def _label(ax, text: str, x: float, y: float, color: str = INK_SECONDARY) -> None:
    """
    Place a small direct label to the right of a point, masking lines behind it.

    Args:
        ax (matplotlib.axes.Axes): Axes to label.
        text (str): Label text.
        x (float): Anchor x in data units.
        y (float): Anchor y in data units.
        color (str): Text colour (an ink token, never a series colour).
    """
    ax.annotate(
        text,
        (x, y),
        xytext=(6, 0),
        textcoords="offset points",
        va="center",
        fontsize=8,
        color=color,
        zorder=4,
        bbox={"boxstyle": "square,pad=0.15", "facecolor": SURFACE, "edgecolor": "none"},
    )


def _legend_above(ax, ncol: int) -> None:
    """
    Put the legend in one row above the plot, clear of the data.

    Args:
        ax (matplotlib.axes.Axes): Axes whose legend to place.
        ncol (int): Number of legend columns.
    """
    ax.legend(
        loc="lower left",
        bbox_to_anchor=(0, 1.02),
        ncol=ncol,
        frameon=False,
        fontsize=8,
        labelcolor=INK_SECONDARY,
        borderaxespad=0,
    )


def _style_axes(ax) -> None:
    """
    Apply the recessive chart chrome: surface fill, hairline axes, muted ticks.

    Args:
        ax (matplotlib.axes.Axes): Axes to style.
    """
    ax.set_facecolor(SURFACE)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(BASELINE)
    ax.tick_params(colors=INK_MUTED, labelcolor=INK_SECONDARY, labelsize=8)
    ax.set_axisbelow(True)


def _save(fig: Figure, path_stem: Path) -> list[Path]:
    """
    Write a figure in every output format.

    Args:
        fig (Figure): Figure to save.
        path_stem (Path): Output path without extension.

    Returns:
        list[Path]: Files written.
    """
    path_stem.parent.mkdir(parents=True, exist_ok=True)
    paths = [path_stem.with_suffix(f".{fmt}") for fmt in FORMATS]
    for path in paths:
        # Tight bbox keeps labels placed outside the axes (end labels, legend) in frame.
        fig.savefig(path, dpi=200, facecolor=SURFACE, bbox_inches="tight", pad_inches=0.1)
    return paths
