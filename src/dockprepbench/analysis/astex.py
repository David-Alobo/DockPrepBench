"""Independent check on the Astex Diverse Set, as specified in docs/astex_validation_plan.md.

H1: the 2 Å box raises top-1 success over the baseline (paired difference with interval,
    exact McNemar, Holm over the four engines).
H2: agreement with the five primary arms has a higher AUROC than the docking score (and, for
    Gnina, than its CNN score); supported if the paired bootstrap interval excludes zero.
H3: the same for the pair chosen on PoseBusters (2 Å box + RDKit conformer).
H4: success when all five arms agree, with a Wilson interval, against the PoseBusters rate;
    descriptive (overlapping intervals read as consistent).
An engine with fewer than 10 incorrect baseline poses has its AUROCs reported but not tested.

Run with `python -m dockprepbench.analysis.astex ASTEX_RESULTS POSEBUSTERS_LEVELS OUT`.
"""

import argparse
from collections.abc import Sequence
from pathlib import Path

import pandas as pd
from rdkit import RDLogger

from dockprepbench.analysis.agreement import bootstrap_aurocs, collect
from dockprepbench.analysis.report import discover
from dockprepbench.analysis.stats import (
    compare,
    holm_adjust,
    paired_difference_interval,
    success_rate,
    successes,
)
from dockprepbench.evaluation.rmsd import SUCCESS_THRESHOLD

BOX_ARM = "box_padding_2"
PAIR = ("box_padding_2", "ligand_rdkit_conformer")  # chosen on PoseBusters before these runs
MIN_NEGATIVES = 10


def h1_table(results: str | Path) -> pd.DataFrame:
    """
    H1: the 2 Å box against the baseline for each engine, Holm-adjusted over the engines.

    Args:
        results (str | Path): Root containing `<engine>/<arm>/` Astex runs.

    Returns:
        pd.DataFrame: One row per engine: rates, difference with 95% paired interval,
            gained, lost, raw and Holm p, and whether H1 is supported.
    """
    rows = []
    for e, r in discover(results).items():
        base, arm = successes(r.baseline), successes(r.arms[BOX_ARM])
        c = compare(base, arm)
        low, high = paired_difference_interval(c.both, len(c.only_a), len(c.only_b), c.n)
        rows.append({
            "engine": e,
            "n": c.n,
            "baseline": base.mean(),
            "box2": arm.mean(),
            "diff_pts": 100 * (arm.mean() - base.mean()),
            "diff_low": 100 * low,
            "diff_high": 100 * high,
            "gained": len(c.only_b),
            "lost": len(c.only_a),
            "p": c.p_value,
        })  # fmt: skip
    table = pd.DataFrame(rows)
    table["p_holm"] = holm_adjust(table["p"].tolist())
    table["supported"] = (table["diff_pts"] > 0) & (table["p_holm"] < 0.05)
    return table


def h2_h3_table(per_complex: pd.DataFrame) -> pd.DataFrame:
    """
    H2 and H3: AUROC of agreement (five arms; the pre-selected pair) against the docking score.

    Args:
        per_complex (pd.DataFrame): `agreement.collect` output for the Astex runs.

    Returns:
        pd.DataFrame: One row per (engine, signal): AUROC with interval, the paired
            difference from the score (and, for Gnina, from the CNN score), the number of
            incorrect baseline poses, and whether the hypothesis is tested and supported.
    """
    rows = []
    for e, g in per_complex.groupby("engine", sort=False):
        g = g.set_index("id")
        negatives = int((~g["success"]).sum())
        signals = pd.DataFrame({
            "score": -g["score"],
            "prep_agree": g["prep_agree"],
            "pair_agree": sum((g[f"d_{a}"] < SUCCESS_THRESHOLD).astype(int) for a in PAIR),
        })  # fmt: skip
        if g["cnn"].notna().any():
            signals["cnn"] = g["cnn"]
        vs_score = bootstrap_aurocs(signals, g["success"], reference="score")
        vs_cnn = (
            bootstrap_aurocs(signals, g["success"], reference="cnn") if "cnn" in signals else None
        )
        for signal, hypothesis in (("prep_agree", "H2"), ("pair_agree", "H3")):
            t = vs_score.loc[signal]
            tested = negatives >= MIN_NEGATIVES and int(g["success"].sum()) >= MIN_NEGATIVES
            row = {
                "engine": e,
                "hypothesis": hypothesis,
                "signal": signal,
                "n": len(g),
                "incorrect": negatives,
                "auroc": t["auroc"],
                "low": t["low"],
                "high": t["high"],
                "score_auroc": vs_score.loc["score", "auroc"],
                "diff_vs_score": t["diff"],
                "diff_low": t["diff_low"],
                "diff_high": t["diff_high"],
                "tested": tested,
                "supported": bool(tested and t["diff_low"] > 0),
            }
            if vs_cnn is not None:
                c = vs_cnn.loc[signal]
                row |= {
                    "cnn_auroc": vs_cnn.loc["cnn", "auroc"],
                    "diff_vs_cnn": c["diff"],
                    "diff_vs_cnn_low": c["diff_low"],
                    "diff_vs_cnn_high": c["diff_high"],
                }
            rows.append(row)
    return pd.DataFrame(rows)


def h4_table(per_complex: pd.DataFrame, posebusters_levels: pd.DataFrame) -> pd.DataFrame:
    """
    H4: success when all five arms agree, on Astex and on PoseBusters.

    Args:
        per_complex (pd.DataFrame): `agreement.collect` output for the Astex runs.
        posebusters_levels (pd.DataFrame): `agreement.level_table` output for PoseBusters.

    Returns:
        pd.DataFrame: One row per engine: unanimous complexes and their success rate with
            Wilson interval on each dataset, and whether the intervals overlap.
    """
    pb = posebusters_levels[posebusters_levels["prep_agree"] == 5].set_index("engine")
    rows = []
    for e, g in per_complex.groupby("engine", sort=False):
        unanimous = g[g["prep_agree"] == 5]
        row = {"engine": e, "n": len(unanimous)}
        if len(unanimous):
            r = success_rate(unanimous["success"])
            row |= {"rate": r.value, "low": r.low, "high": r.high}
        if e in pb.index:
            p = pb.loc[e]
            row |= {
                "pb_n": int(p["n"]),
                "pb_rate": p["rate"],
                "pb_low": p["low"],
                "pb_high": p["high"],
            }
            if len(unanimous):
                row["consistent"] = bool(r.low <= p["high"] and p["low"] <= r.high)
        rows.append(row)
    return pd.DataFrame(rows)


def main(argv: Sequence[str] | None = None) -> None:
    """
    Run H1-H4 on the Astex runs and write the tables as CSV.

    Args:
        argv (Sequence[str] | None): Command-line arguments; sys.argv if None.
    """
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("results", type=Path, help="Astex root with <engine>/<arm>/ runs")
    parser.add_argument("levels", type=Path, help="PoseBusters agreement_levels.csv")
    parser.add_argument("out", type=Path, help="directory for the CSV tables")
    args = parser.parse_args(argv)
    RDLogger.DisableLog("rdApp.*")
    args.out.mkdir(parents=True, exist_ok=True)
    per_complex = collect(args.results)
    baseline_rates = {
        e: success_rate(successes(r.baseline)) for e, r in discover(args.results).items()
    }
    tables = {
        "astex_baseline": pd.DataFrame(
            [
                {"engine": e, "n": r.total, "rate": r.value, "low": r.low, "high": r.high}
                for e, r in baseline_rates.items()
            ]  # fmt: skip
        ),
        "astex_h1": h1_table(args.results),
        "astex_h2_h3": h2_h3_table(per_complex),
        "astex_h4": h4_table(per_complex, pd.read_csv(args.levels)),
    }
    for name, table in tables.items():
        table.to_csv(args.out / f"{name}.csv", index=False)
        print(f"\n== {name}\n{table.round(3).to_string(index=False)}")


if __name__ == "__main__":
    main()
