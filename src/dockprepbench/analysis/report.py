"""Benchmark report: every table in the paper, regenerated from the raw run results.

Expects runs laid out by `dockprepbench run`: `<results>/<engine>/<arm>/results.csv`.
A run named `baseline_seed<N>` is a seed replicate of that engine's baseline.
"""

import itertools
import re
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from dockprepbench.analysis.stats import compare, holm_adjust, load_run, success_rate, successes

BASELINE = "baseline"
SEED_RUN = re.compile(rf"^{BASELINE}_seed(\d+)$")
SITE_SHIFT = 4.0  # Å; heuristic split between wrong site and wrong binding mode
METRICS = {"top1": False, "pb_valid": True}  # metric name -> require PB-valid
CLASSES = ["valid", "invalid", "scoring", "site", "mode", "failed"]


@dataclass(frozen=True)
class EngineRuns:
    """One engine's runs: the baseline, its seed replicates and the preparation arms."""

    engine: str
    baseline: pd.DataFrame
    seeds: dict[str, pd.DataFrame]  # replicate name -> results, excluding the baseline itself
    arms: dict[str, pd.DataFrame]


def discover(results: str | Path) -> dict[str, EngineRuns]:
    """
    Find every finished run under a results directory, grouped by engine.

    Engines without a finished baseline are skipped, since every comparison is to it.

    Args:
        results (str | Path): Root containing `<engine>/<arm>/results.csv`.

    Returns:
        dict[str, EngineRuns]: Runs per engine, in sorted order.
    """
    found = {}
    for engine_dir in sorted(p for p in Path(results).iterdir() if p.is_dir()):
        runs = {
            d.name: load_run(d / "results.csv")
            for d in sorted(engine_dir.iterdir())
            if (d / "results.csv").is_file()
        }
        if BASELINE not in runs:
            continue
        seeds = {k: v for k, v in runs.items() if SEED_RUN.match(k)}
        arms = {k: v for k, v in runs.items() if k != BASELINE and k not in seeds}
        found[engine_dir.name] = EngineRuns(engine_dir.name, runs[BASELINE], seeds, arms)
    return found


def restrict(runs: dict[str, EngineRuns], ids: Iterable[str]) -> dict[str, EngineRuns]:
    """
    Keep only the given complexes in every run.

    Args:
        runs (dict[str, EngineRuns]): Output of `discover`.
        ids (Iterable[str]): Complex IDs to keep; IDs absent from a run are ignored.

    Returns:
        dict[str, EngineRuns]: The same runs, each limited to `ids`.
    """
    keep = list(dict.fromkeys(ids))

    def cut(run: pd.DataFrame) -> pd.DataFrame:
        return run.loc[run.index.intersection(keep)]

    return {
        name: EngineRuns(
            e.engine,
            cut(e.baseline),
            {k: cut(v) for k, v in e.seeds.items()},
            {k: cut(v) for k, v in e.arms.items()},
        )
        for name, e in runs.items()
    }


def read_ids(path: str | Path) -> list[str]:
    """
    Read a subset file: one complex ID per line; blank lines and `#` comments are ignored.

    Args:
        path (str | Path): Subset file, e.g. configs/subsets/posebusters_308.txt.

    Returns:
        list[str]: Complex IDs in file order.
    """
    lines = (line.split("#")[0].strip() for line in Path(path).read_text().splitlines())
    return [line for line in lines if line]


def outcome_class(row: pd.Series, site_shift: float = SITE_SHIFT) -> str:
    """
    Classify one complex's top-1 outcome (see the pose-failure-triage conventions).

    Args:
        row (pd.Series): One row of a results table.
        site_shift (float): Centroid shift in Å above which a miss is a wrong site.

    Returns:
        str: One of CLASSES.
    """
    if row["error"]:
        return "failed"
    if row["top1_rmsd"] < 2:
        return "valid" if row["top1_pb_valid"] else "invalid"
    if row["best_rmsd"] < 2:
        return "scoring"
    return "site" if row["top1_centroid_shift"] > site_shift else "mode"


def summary_table(runs: dict[str, EngineRuns]) -> pd.DataFrame:
    """
    Success rates with 95% Wilson intervals for every run.

    Args:
        runs (dict[str, EngineRuns]): Output of `discover`.

    Returns:
        pd.DataFrame: One row per engine and run.
    """
    rows = []
    for e in runs.values():
        for name, run in {BASELINE: e.baseline, **e.arms, **e.seeds}.items():
            row = {
                "engine": e.engine,
                "run": name,
                "n": len(run),
                "failed": int((run["error"] != "").sum()),
            }
            for metric, pb in METRICS.items():
                rate = success_rate(successes(run, pb_valid=pb))
                row |= {metric: rate.value, f"{metric}_low": rate.low, f"{metric}_high": rate.high}
            rows.append(row)
    return pd.DataFrame(rows)


def comparison_table(runs: dict[str, EngineRuns]) -> pd.DataFrame:
    """
    Paired comparison of each arm with its engine's baseline, Holm-adjusted per engine and metric.

    Args:
        runs (dict[str, EngineRuns]): Output of `discover`.

    Returns:
        pd.DataFrame: One row per engine, metric and arm.
    """
    rows = []
    for e in runs.values():
        for metric, pb in METRICS.items():
            base = successes(e.baseline, pb_valid=pb)
            family = []
            for arm, run in e.arms.items():
                hits = successes(run, pb_valid=pb)
                c = compare(base, hits)
                shared = base.index.intersection(hits.index)
                family.append({
                    "engine": e.engine,
                    "metric": metric,
                    "arm": arm,
                    "n": c.n,
                    "baseline": base[shared].mean(),
                    "arm_rate": hits[shared].mean(),
                    "delta_pts": 100 * (hits[shared].mean() - base[shared].mean()),
                    "gained": len(c.only_b),
                    "lost": len(c.only_a),
                    "p": c.p_value,
                })  # fmt: skip
            for row, p_holm in zip(family, holm_adjust([r["p"] for r in family]), strict=True):
                rows.append(row | {"p_holm": p_holm, "family_size": len(family)})
    return pd.DataFrame(rows)


def noise_table(runs: dict[str, EngineRuns]) -> pd.DataFrame:
    """
    Seed-to-seed variation of each engine's baseline: rates and per-complex outcome changes.

    Args:
        runs (dict[str, EngineRuns]): Output of `discover`.

    Returns:
        pd.DataFrame: One row per engine with at least one seed replicate.
    """
    rows = []
    for e in runs.values():
        if not e.seeds:
            continue
        replicates = {BASELINE: e.baseline, **e.seeds}
        hits = {k: successes(v) for k, v in replicates.items()}
        rates = pd.Series({k: 100 * h.mean() for k, h in hits.items()})
        flips = [int((hits[a] != hits[b]).sum()) for a, b in itertools.combinations(hits, 2)]
        table = pd.DataFrame(hits)
        always, never = int(table.all(axis=1).sum()), int((~table).all(axis=1).sum())
        rows.append({
            "engine": e.engine,
            "replicates": len(hits),
            "top1_mean": rates.mean(),
            "top1_sd": rates.std(),
            "top1_range": rates.max() - rates.min(),
            "flips_mean": sum(flips) / len(flips),
            "flips_min": min(flips),
            "flips_max": max(flips),
            "always_succeed": always,
            "always_fail": never,
            "seed_sensitive": len(table) - always - never,
        })  # fmt: skip
    return pd.DataFrame(rows)


def flips_table(runs: dict[str, EngineRuns], noise: pd.DataFrame) -> pd.DataFrame:
    """
    How many complexes each arm flips relative to the baseline, against seed-to-seed flips.

    Args:
        runs (dict[str, EngineRuns]): Output of `discover`.
        noise (pd.DataFrame): Output of `noise_table`.

    Returns:
        pd.DataFrame: One row per engine and arm; `vs_noise` is NaN without seed replicates.
    """
    seed_flips = dict(zip(noise.get("engine", []), noise.get("flips_mean", []), strict=True))
    rows = []
    for e in runs.values():
        base = successes(e.baseline)
        for arm, run in e.arms.items():
            hits = successes(run)
            shared = base.index.intersection(hits.index)
            flipped = int((base[shared] != hits[shared]).sum())
            noise_level = seed_flips.get(e.engine, float("nan"))
            rows.append({
                "engine": e.engine,
                "arm": arm,
                "delta_pts": 100 * (hits[shared].mean() - base[shared].mean()),
                "flipped": flipped,
                "seed_flips_mean": noise_level,
                "vs_noise": flipped / noise_level if noise_level else float("nan"),
            })  # fmt: skip
    return pd.DataFrame(rows)


def transition_table(runs: dict[str, EngineRuns], site_shift: float = SITE_SHIFT) -> pd.DataFrame:
    """
    Counts of complexes moving between outcome classes, baseline -> arm.

    Args:
        runs (dict[str, EngineRuns]): Output of `discover`.
        site_shift (float): Centroid shift threshold for the wrong-site class.

    Returns:
        pd.DataFrame: One row per engine, arm and (from, to) pair with a non-zero count.
    """
    rows = []
    for e in runs.values():
        base = e.baseline.apply(outcome_class, axis=1, site_shift=site_shift)
        for arm, run in e.arms.items():
            arm_classes = run.apply(outcome_class, axis=1, site_shift=site_shift)
            shared = base.index.intersection(arm_classes.index)
            pairs = pd.crosstab(base[shared], arm_classes[shared]).stack()
            for (src, dst), count in pairs.items():
                if count:
                    rows.append(
                        {"engine": e.engine, "arm": arm, "from": src, "to": dst, "n": int(count)}
                    )
    return pd.DataFrame(rows)


def write_report(
    results: str | Path,
    out: str | Path,
    dataset: str | Path | None = None,
    subset: Iterable[str] | None = None,
) -> list[Path]:
    """
    Regenerate all report tables from raw results and write them as CSV files.

    Args:
        results (str | Path): Root containing `<engine>/<arm>/results.csv`.
        out (str | Path): Directory for the tables.
        dataset (str | Path | None): Extracted PoseBusters set; when given, also writes the
            ligand properties and the sensitivity of each arm's flips to them.
        subset (Iterable[str] | None): Complex IDs to restrict every table to.

    Returns:
        list[Path]: Files written.

    Raises:
        ValueError: If no engine has a finished baseline.
    """
    runs = discover(results)
    if not runs:
        raise ValueError(f"No finished baseline runs under {results}")
    if subset is not None:
        runs = restrict(runs, subset)
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    noise = noise_table(runs)
    tables = {
        "summary": summary_table(runs),
        "comparisons": comparison_table(runs),
        "noise": noise,
        "flips": flips_table(runs, noise),
        "transitions": transition_table(runs),
    }
    if dataset is not None:
        # Imported here: the sensitivity module builds on this one.
        from dockprepbench.analysis.sensitivity import ligand_properties, sensitivity_table
        from dockprepbench.datasets.posebusters import load_complexes

        properties = ligand_properties(load_complexes(dataset))
        tables["ligand_properties"] = properties.reset_index()
        tables["sensitivity"] = sensitivity_table(runs, properties)
    written = []
    for name, table in tables.items():
        path = out / f"{name}.csv"
        table.to_csv(path, index=False, float_format="%.6g")
        written.append(path)
    # LaTeX fragments for the supplementary information, from the CSVs just written.
    from dockprepbench.analysis.latex import write_latex

    return written + write_latex(out)
