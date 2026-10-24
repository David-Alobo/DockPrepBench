"""LaTeX table bodies for the supplementary information, generated from the report tables.

Each function returns the rows of a `tabular` or `longtable` (header included), so the
document keeps control of captions and placement while every number comes from the
report. Engines and arms appear in the same order and with the same names as in the
figures.
"""

import math
from pathlib import Path

import pandas as pd

from dockprepbench.analysis.figures import ARM_LABELS, ENGINE_LABELS

METRIC_LABELS = {"top1": "top-1", "pb_valid": "PB-valid"}
PROPERTY_LABELS = {
    "heavy_atoms": "Heavy atoms",
    "rotatable_bonds": "Rotatable bonds",
    "abs_charge": "Absolute charge",
    "charge_changed": "Charge changed at pH 7.4",
    "cofactor": "Cofactor present",
    "seed_sensitive": "Seed-sensitive",
}


def _arm(key: str) -> str:
    """
    Name an arm for a table, with LaTeX-safe units.

    Args:
        key (str): Arm key, e.g. box_padding_2.

    Returns:
        str: Display name.
    """
    return ARM_LABELS.get(key, key.replace("_", r"\_")).replace("Å", r"\AA{}")


def _engines(frame: pd.DataFrame) -> list[str]:
    """
    Order the engines in a table: known engines first, as in the figures.

    Args:
        frame (pd.DataFrame): Table with an `engine` column.

    Returns:
        list[str]: Engine keys in display order.
    """
    found = set(frame["engine"])
    return [e for e in ENGINE_LABELS if e in found] + sorted(found - set(ENGINE_LABELS))


def _pct(x: float) -> str:
    """
    Format a fraction as a percentage with one decimal.

    Args:
        x (float): Fraction.

    Returns:
        str: E.g. "46.3".
    """
    return f"{100 * x:.1f}"


def _p(p: float) -> str:
    """
    Format a p-value: two significant figures, scientific below 0.001.

    Args:
        p (float): P-value.

    Returns:
        str: LaTeX math for the value.
    """
    if not math.isfinite(p):
        return "--"
    if p >= 0.001:
        return f"{p:#.2g}" if p < 0.995 else "1.0"  # keeps trailing zeros: 0.20
    mantissa, exponent = f"{p:.1e}".split("e")
    return rf"${mantissa}\times10^{{{int(exponent)}}}$"


def _signed(x: float) -> str:
    """
    Format a change in points with an explicit sign (typographic minus).

    Args:
        x (float): Change in percentage points.

    Returns:
        str: E.g. "+6.8" or "$-1.9$".
    """
    return f"+{x:.1f}" if x >= 0 else f"${x:.1f}$"


def arm_table(summary: pd.DataFrame, comparisons: pd.DataFrame, metric: str) -> str:
    """
    Success per engine and arm with intervals, and each arm's paired comparison.

    Args:
        summary (pd.DataFrame): The report's summary table.
        comparisons (pd.DataFrame): The report's comparisons table.
        metric (str): "top1" or "pb_valid".

    Returns:
        str: Header and rows for a 7-column tabular (l l r r r r r).
    """
    rate = summary.set_index(["engine", "run"])
    comp = comparisons[comparisons["metric"] == metric].set_index(["engine", "arm"])
    lines = [
        r"\toprule",
        rf"Engine & Arm & {METRIC_LABELS[metric].capitalize()} (\%) [95\% CI] & $\Delta$ (pts) "
        r"& Gained & Lost & Holm $p$ \\",
        r"\midrule",
    ]
    for i, engine in enumerate(_engines(summary)):
        if i:
            lines.append(r"\addlinespace")
        arms = [a for a in ARM_LABELS if (engine, a) in rate.index]
        for j, arm in enumerate(arms):
            r = rate.loc[(engine, arm)]
            cell = f"{_pct(r[metric])} [{_pct(r[f'{metric}_low'])}, {_pct(r[f'{metric}_high'])}]"
            name = ENGINE_LABELS.get(engine, engine) if j == 0 else ""
            if arm == "baseline" or (engine, arm) not in comp.index:
                lines.append(rf"{name} & {_arm(arm)} & {cell} & -- & -- & -- & -- \\")
                continue
            c = comp.loc[(engine, arm)]
            lines.append(
                rf"{name} & {_arm(arm)} & {cell} & {_signed(c['delta_pts'])} & {int(c['gained'])} "
                rf"& {int(c['lost'])} & {_p(c['p_holm'])} \\"
            )
    lines.append(r"\bottomrule")
    return "\n".join(lines) + "\n"


def flips_table(flips: pd.DataFrame, noise: pd.DataFrame) -> str:
    """
    Complexes each arm flips relative to its baseline, against seed-to-seed flips.

    Args:
        flips (pd.DataFrame): The report's flips table.
        noise (pd.DataFrame): The report's noise table.

    Returns:
        str: Header and rows for a 5-column tabular (l l r r r).
    """
    seeds = noise.set_index("engine") if len(noise) else pd.DataFrame()
    lines = [
        r"\toprule",
        r"Engine & Arm & $\Delta$ (pts) & Flipped & Ratio to seed flips \\",
        r"\midrule",
    ]
    table = flips.set_index(["engine", "arm"])
    for i, engine in enumerate(_engines(flips)):
        if i:
            lines.append(r"\addlinespace")
        arms = [a for a in ARM_LABELS if (engine, a) in table.index]
        for j, arm in enumerate(arms):
            r = table.loc[(engine, arm)]
            name = ENGINE_LABELS.get(engine, engine) if j == 0 else ""
            ratio = "--" if engine not in seeds.index else f"{r['vs_noise']:.1f}"
            lines.append(
                rf"{name} & {_arm(arm)} & {_signed(r['delta_pts'])} & {int(r['flipped'])} "
                rf"& {ratio} \\"
            )
    lines.append(r"\bottomrule")
    return "\n".join(lines) + "\n"


def sensitivity_table(sensitivity: pd.DataFrame) -> str:
    """
    Flipped against unflipped complexes on each ligand property, per engine and arm.

    Args:
        sensitivity (pd.DataFrame): The report's sensitivity table.

    Returns:
        str: Header and rows for a 7-column longtable (l l l r r r r), header repeated
            on each page.
    """
    header = [
        r"\toprule",
        r"Engine & Arm & Property & Flipped & Unflipped & Effect & Holm $p$ \\",
        r"\midrule",
    ]
    lines = [*header, r"\endfirsthead", *header, r"\endhead", r"\bottomrule", r"\endfoot"]
    table = sensitivity.set_index(["engine", "arm", "property"])
    for engine in _engines(sensitivity):
        for arm in [a for a in ARM_LABELS if (engine, a) in table.droplevel(2).index]:
            rows = table.loc[(engine, arm)]
            for k, (prop, r) in enumerate(rows.iterrows()):
                name = ENGINE_LABELS.get(engine, engine) if k == 0 else ""
                arm_name = _arm(arm) if k == 0 else ""
                if r["effect_kind"] == "odds ratio":
                    values = f"{_pct(r['flipped_value'])}\\% & {_pct(r['unflipped_value'])}\\%"
                    effect = "$\\infty$" if math.isinf(r["effect"]) else f"OR {r['effect']:.2g}"
                else:
                    values = f"{r['flipped_value']:g} & {r['unflipped_value']:g}"
                    effect = f"$r = {r['effect']:.2f}$"
                lines.append(
                    rf"{name} & {arm_name} & {PROPERTY_LABELS.get(prop, prop)} & {values} "
                    rf"& {effect} & {_p(r['p_holm'])} \\"
                )
            lines.append(r"\addlinespace")
    return "\n".join(lines) + "\n"


def write_latex(tables: Path) -> list[Path]:
    """
    Write the supplementary tables as LaTeX fragments next to the report's CSV files.

    Args:
        tables (Path): Directory written by `report.write_report`.

    Returns:
        list[Path]: Fragments written (`*.tex`); tables whose CSV is missing are skipped.
    """
    read = lambda name: pd.read_csv(tables / f"{name}.csv")  # noqa: E731

    def tabular(spec: str, body: str) -> str:
        # The whole environment goes in the file: LaTeX cannot \input rows into a tabular.
        return f"\\begin{{tabular}}{{{spec}}}\n{body}\\end{{tabular}}\n"

    def longtable(spec: str, caption: str, label: str, body: str) -> str:
        # The caption text is a macro the including document defines, so captions stay there.
        return (
            f"\\begin{{longtable}}{{{spec}}}\n\\caption{{{caption}}}\\label{{{label}}}\\\\\n"
            f"{body}\\end{{longtable}}\n"
        )

    jobs = {
        "top1_table": lambda: tabular(
            "llrrrrr", arm_table(read("summary"), read("comparisons"), "top1")
        ),
        "pb_valid_table": lambda: tabular(
            "llrrrrr", arm_table(read("summary"), read("comparisons"), "pb_valid")
        ),
        "flips_table": lambda: tabular("llrrr", flips_table(read("flips"), read("noise"))),
        "sensitivity_table": lambda: longtable(
            "lllrrrr",
            r"\SensitivityCaption",
            "tab:si-sensitivity",
            sensitivity_table(read("sensitivity")),
        ),
    }
    written = []
    for name, make in jobs.items():
        try:
            body = make()
        except FileNotFoundError:
            continue
        path = tables / f"{name}.tex"
        path.write_text(body, encoding="utf-8")
        written.append(path)
    return written
