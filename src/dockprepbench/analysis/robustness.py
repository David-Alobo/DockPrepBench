"""Robustness checks requested in review: the wrong-site threshold and the removed cofactors.

1. The failure-mode split between a wrong site and a wrong binding mode uses a 4 Å centroid
   shift. Here the wrong-site counts are recomputed at other thresholds, to show whether the
   box-size conclusion depends on that choice.
2. Every heterogen other than additives is removed from the receptor before docking, while
   the PoseBusters checks use the original protein with its cofactors. Here the complexes
   with cofactors are identified, the PoseBusters failures caused by cofactors counted, and
   the main comparisons repeated on the complexes without cofactors.

Run with `python -m dockprepbench.analysis.robustness RESULTS DATASET OUT`.
"""

import argparse
from collections.abc import Sequence
from pathlib import Path

import pandas as pd

from dockprepbench.analysis.report import (
    EngineRuns,
    comparison_table,
    discover,
    outcome_class,
    restrict,
    summary_table,
)
from dockprepbench.analysis.sensitivity import cofactors
from dockprepbench.datasets.posebusters import load_complexes

THRESHOLDS = (2.0, 3.0, 4.0, 5.0, 6.0)
COFACTOR_CHECKS = (
    "minimum_distance_to_organic_cofactors",
    "volume_overlap_with_organic_cofactors",
    "minimum_distance_to_inorganic_cofactors",
    "volume_overlap_with_inorganic_cofactors",
)


def threshold_table(
    runs: dict[str, EngineRuns], thresholds: Sequence[float] = THRESHOLDS
) -> pd.DataFrame:
    """
    Top-1 wrong-site and wrong-mode counts per engine and run at several centroid thresholds.

    Args:
        runs (dict[str, EngineRuns]): Output of `discover`.
        thresholds (Sequence[float]): Centroid shifts in Å above which a miss is a wrong site.

    Returns:
        pd.DataFrame: One row per (engine, run, threshold): counts of wrong-site and
            wrong-mode top-1 poses.
    """
    rows = []
    for e in runs.values():
        for name, run in {"baseline": e.baseline, **e.arms}.items():
            for t in thresholds:
                classes = run.apply(outcome_class, axis=1, site_shift=t)
                rows.append({
                    "engine": e.engine,
                    "run": name,
                    "threshold": t,
                    "site": int((classes == "site").sum()),
                    "mode": int((classes == "mode").sum()),
                })  # fmt: skip
    return pd.DataFrame(rows)


def cofactor_ids(dataset: str | Path) -> dict[str, list[str]]:
    """
    Cofactors found in each complex's protein file.

    Args:
        dataset (str | Path): Extracted PoseBusters root (as returned by `download`).

    Returns:
        dict[str, list[str]]: Complex ID to its cofactor residue names (empty if none).
    """
    return {c.id: cofactors(c) for c in load_complexes(dataset)}


def cofactor_failure_table(runs: dict[str, EngineRuns], has_cofactor: pd.Series) -> pd.DataFrame:
    """
    How often a correct top-1 pose fails PoseBusters only because of a removed cofactor.

    Args:
        runs (dict[str, EngineRuns]): Output of `discover`.
        has_cofactor (pd.Series): Boolean per complex ID.

    Returns:
        pd.DataFrame: One row per (engine, run): complexes with a cofactor, top-1 successes
            among them, and successes failing PoseBusters on a cofactor check alone.
    """
    rows = []
    for e in runs.values():
        for name, run in {"baseline": e.baseline, **e.arms}.items():
            hit = run["top1_rmsd"] < 2
            checks = run["top1_failed_checks"].fillna("").str.split(";")
            cofactor_only = checks.map(lambda c: bool(c[0]) and set(c) <= set(COFACTOR_CHECKS))
            with_cof = has_cofactor.reindex(run.index, fill_value=False)
            rows.append({
                "engine": e.engine,
                "run": name,
                "with_cofactor": int(with_cof.sum()),
                "top1_success_with_cofactor": int((hit & with_cof).sum()),
                "failed_on_cofactor_only": int((hit & cofactor_only).sum()),
            })  # fmt: skip
    return pd.DataFrame(rows)


def main(argv: Sequence[str] | None = None) -> None:
    """
    Run both robustness checks and write their tables as CSV.

    Args:
        argv (Sequence[str] | None): Command-line arguments; sys.argv if None.
    """
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("results", type=Path, help="root with <engine>/<arm>/ runs")
    parser.add_argument("dataset", type=Path, help="extracted PoseBusters root")
    parser.add_argument("out", type=Path, help="directory for the CSV tables")
    args = parser.parse_args(argv)
    args.out.mkdir(parents=True, exist_ok=True)
    runs = discover(args.results)
    found = cofactor_ids(args.dataset)
    has_cofactor = pd.Series({cid: bool(c) for cid, c in found.items()})
    without = [cid for cid, c in found.items() if not c]
    subset = restrict(runs, without)
    tables = {
        "robust_threshold": threshold_table(runs),
        "robust_cofactor_list": pd.DataFrame(
            [{"id": cid, "cofactors": " ".join(c)} for cid, c in found.items() if c]
        ),
        "robust_cofactor_failures": cofactor_failure_table(runs, has_cofactor),
        "robust_no_cofactor_summary": summary_table(subset),
        "robust_no_cofactor_comparisons": comparison_table(subset),
    }
    for name, table in tables.items():
        table.to_csv(args.out / f"{name}.csv", index=False)
    print(f"{int(has_cofactor.sum())} of {len(has_cofactor)} complexes have a cofactor")


if __name__ == "__main__":
    main()
