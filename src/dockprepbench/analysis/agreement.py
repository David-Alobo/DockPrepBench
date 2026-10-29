"""Pose agreement across preparation arms as a confidence signal, and the consensus pose.

For each engine and complex, the baseline top-1 pose is compared with the top-1 pose of each
primary arm (symmetry-corrected heavy-atom RMSD, no superposition). The number of arms whose
pose lies within 2 Å of the baseline pose ("preparation agreement", 0-5) is tested as a
predictor of whether the baseline pose is correct, against the docking score, Gnina's CNN
score, agreement across engines and, for Vina, agreement across seed replicates (a control:
same preparation, different random seed).

The consensus pose is the top-1 pose, among the baseline and the arms, that most of the other
candidates agree with; it is compared with the baseline pose complex by complex.

Run with `python -m dockprepbench.analysis.agreement RESULTS OUT`.
"""

import argparse
import itertools
from collections.abc import Sequence
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger

from dockprepbench.analysis.report import BASELINE, PRIMARY_ARMS, discover
from dockprepbench.analysis.stats import compare, holm_adjust, success_rate, successes
from dockprepbench.evaluation.rmsd import SUCCESS_THRESHOLD, pose_rmsd

N_BOOTSTRAP = 2000
SEED = 0


def top_pose(run_dir: Path, cid: str) -> Chem.Mol | None:
    """
    Read the top-ranked pose a run saved for one complex.

    Args:
        run_dir (Path): Run directory containing `poses/<id>.sdf`, poses in rank order.
        cid (str): Complex ID.

    Returns:
        Chem.Mol | None: The first pose, or None if the run saved no poses for the complex.
    """
    path = run_dir / "poses" / f"{cid}.sdf"
    if not path.is_file():
        return None
    supplier = Chem.SDMolSupplier(str(path), removeHs=False, sanitize=False)
    return next((mol for mol in supplier if mol is not None), None)


def pair_rmsd(a: Chem.Mol | None, b: Chem.Mol | None) -> float:
    """
    RMSD between two poses of the same ligand, or NaN if either is missing or they don't match.

    Args:
        a (Chem.Mol | None): First pose.
        b (Chem.Mol | None): Second pose.

    Returns:
        float: Symmetry-corrected heavy-atom RMSD in Å, or NaN.
    """
    if a is None or b is None:
        return float("nan")
    try:
        return pose_rmsd(a, b)
    except (RuntimeError, ValueError):
        return float("nan")


def agreement(distances: pd.DataFrame, cutoff: float = SUCCESS_THRESHOLD) -> pd.Series:
    """
    Count, per row, the distances below a cutoff. Missing distances count as disagreement.

    Args:
        distances (pd.DataFrame): RMSDs from a reference pose, one column per other pose.
        cutoff (float): Agreement threshold in Å.

    Returns:
        pd.Series: Integer count per row.
    """
    return (distances < cutoff).sum(axis=1).astype(int)


def auroc(signal: Sequence[float], labels: Sequence[bool]) -> float:
    """
    Area under the ROC curve: the chance that a random positive outranks a random negative.

    Ties count one half. Rows with a NaN signal are dropped.

    Args:
        signal (Sequence[float]): Higher means more confident.
        labels (Sequence[bool]): True for positives.

    Returns:
        float: AUROC in [0, 1], or NaN if either class is empty.
    """
    s = np.asarray(signal, dtype=float)
    y = np.asarray(labels, dtype=bool)
    keep = ~np.isnan(s)
    s, y = s[keep], y[keep]
    n_pos, n_neg = int(y.sum()), int((~y).sum())
    if n_pos == 0 or n_neg == 0:
        return float("nan")
    ranks = pd.Series(s).rank(method="average").to_numpy()
    return float((ranks[y].sum() - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg))


def bootstrap_aurocs(
    signals: pd.DataFrame, labels: pd.Series, reference: str, n: int = N_BOOTSTRAP, seed: int = SEED
) -> pd.DataFrame:
    """
    AUROC of each signal with a 95% percentile interval, and its difference from a reference.

    Complexes are resampled with replacement; all signals share the resamples, so the
    differences are paired.

    Args:
        signals (pd.DataFrame): One column per signal (higher = more confident), indexed by complex.
        labels (pd.Series): Success flag per complex, same index.
        reference (str): Column the differences are taken against.
        n (int): Number of bootstrap resamples.
        seed (int): Random seed.

    Returns:
        pd.DataFrame: One row per signal: auroc, low, high, diff, diff_low, diff_high.
    """
    rng = np.random.default_rng(seed)
    y = labels.to_numpy(dtype=bool)
    values = {c: signals[c].to_numpy(dtype=float) for c in signals}
    point = {c: auroc(v, y) for c, v in values.items()}
    boot = {c: [] for c in values}
    for _ in range(n):
        idx = rng.integers(0, len(y), len(y))
        for c, v in values.items():
            boot[c].append(auroc(v[idx], y[idx]))
    rows = []
    for c in values:
        b = np.array(boot[c])
        d = b - np.array(boot[reference])
        rows.append({
            "signal": c,
            "auroc": point[c],
            "low": np.nanpercentile(b, 2.5),
            "high": np.nanpercentile(b, 97.5),
            "diff": point[c] - point[reference],
            "diff_low": np.nanpercentile(d, 2.5),
            "diff_high": np.nanpercentile(d, 97.5),
        })  # fmt: skip
    return pd.DataFrame(rows).set_index("signal")


def consensus_index(poses: Sequence[Chem.Mol | None], cutoff: float = SUCCESS_THRESHOLD) -> int:
    """
    Pick the pose that most of the other poses agree with.

    Ties go to the earliest pose, so with the baseline first a tie keeps the baseline pose.
    Missing poses are never picked unless every pose is missing.

    Args:
        poses (Sequence[Chem.Mol | None]): Candidate poses, the baseline first.
        cutoff (float): Agreement threshold in Å.

    Returns:
        int: Index of the chosen pose.
    """
    votes = [0 if p is None else 1 for p in poses]  # a present pose agrees with itself
    for i, j in itertools.combinations(range(len(poses)), 2):
        if pair_rmsd(poses[i], poses[j]) < cutoff:
            votes[i] += 1
            votes[j] += 1
    return int(np.argmax(votes))


def collect(results: str | Path, engines: Sequence[str] | None = None) -> pd.DataFrame:
    """
    Compute every per-complex quantity the analysis needs, for each engine.

    Args:
        results (str | Path): Root containing `<engine>/<arm>/` run directories.
        engines (Sequence[str] | None): Engines to include; all with a baseline if None.

    Returns:
        pd.DataFrame: One row per (engine, complex) with baseline success, PB-valid success,
            score, CNN score, preparation, engine and seed agreement, and the consensus
            pose's arm, success and PB-valid success.
    """
    root = Path(results)
    runs = discover(root)
    engines = [e for e in (engines or runs) if e in runs]
    candidates = [BASELINE, *PRIMARY_ARMS]
    poses = {
        (e, a): {cid: top_pose(root / e / a, cid) for cid in runs[e].baseline.index}
        for e in engines
        for a in candidates + list(runs[e].seeds)
    }
    rows = []
    for e in engines:
        r = runs[e]
        frames = {BASELINE: r.baseline, **r.arms}
        hit = {a: successes(f) for a, f in frames.items()}
        valid = {a: successes(f, pb_valid=True) for a, f in frames.items()}
        for cid in r.baseline.index:
            base = poses[(e, BASELINE)][cid]
            arm_d = {a: pair_rmsd(base, poses[(e, a)][cid]) for a in PRIMARY_ARMS}
            eng_d = {o: pair_rmsd(base, poses[(o, BASELINE)][cid]) for o in engines if o != e}
            seed_d = {s: pair_rmsd(base, poses[(e, s)][cid]) for s in r.seeds}
            seed_agree = sum(d < SUCCESS_THRESHOLD for d in seed_d.values()) if seed_d else np.nan
            chosen = candidates[consensus_index([poses[(e, a)][cid] for a in candidates])]
            cnn = float("nan")
            if base is not None and base.HasProp("CNNscore"):
                cnn = float(base.GetProp("CNNscore"))
            rows.append({
                "engine": e,
                "id": cid,
                "success": bool(hit[BASELINE][cid]),
                "pb_success": bool(valid[BASELINE][cid]),
                "score": r.baseline.at[cid, "top1_score"],
                "cnn": cnn,
                "prep_agree": int(sum(d < SUCCESS_THRESHOLD for d in arm_d.values())),
                "engine_agree": int(sum(d < SUCCESS_THRESHOLD for d in eng_d.values())),
                "seed_agree": seed_agree,
                **{f"d_{a}": d for a, d in arm_d.items()},
                "consensus_arm": chosen,
                "consensus_success": bool(hit[chosen].get(cid, False)),
                "consensus_pb_success": bool(valid[chosen].get(cid, False)),
            })  # fmt: skip
    return pd.DataFrame(rows)


def signal_table(per_complex: pd.DataFrame) -> pd.DataFrame:
    """
    AUROC of each confidence signal for baseline top-1 success, with bootstrap intervals.

    Differences are against the docking score and against preparation agreement. Lower docking
    scores are better, so the score enters negated.

    Args:
        per_complex (pd.DataFrame): Output of `collect`.

    Returns:
        pd.DataFrame: One row per (engine, signal).
    """
    tables = []
    for e, g in per_complex.groupby("engine", sort=False):
        g = g.set_index("id")
        signals = pd.DataFrame({
            "score": -g["score"],
            "prep_agree": g["prep_agree"],
            "engine_agree": g["engine_agree"],
            "prep_plus_engine": g["prep_agree"] + g["engine_agree"],
        })  # fmt: skip
        if g["cnn"].notna().any():
            signals["cnn"] = g["cnn"]
        if g["seed_agree"].notna().any():
            signals["seed_agree"] = g["seed_agree"]
        t = bootstrap_aurocs(signals, g["success"], reference="score")
        # Same seed, so the same resamples: differences from preparation agreement are paired too.
        vs_prep = bootstrap_aurocs(signals, g["success"], reference="prep_agree")
        t = t.join(vs_prep[["diff", "diff_low", "diff_high"]].add_suffix("_vs_prep"))
        tables.append(t.assign(engine=e, n=len(g)).reset_index())
    return pd.concat(tables, ignore_index=True)


def level_table(per_complex: pd.DataFrame) -> pd.DataFrame:
    """
    Baseline top-1 success at each level of preparation agreement, with Wilson intervals.

    Args:
        per_complex (pd.DataFrame): Output of `collect`.

    Returns:
        pd.DataFrame: One row per (engine, agreement level): n, success rate and interval.
    """
    rows = []
    for (e, level), g in per_complex.groupby(["engine", "prep_agree"], sort=False):
        rate = success_rate(g["success"])
        rows.append({"engine": e, "prep_agree": level, "n": rate.total, "rate": rate.value,
                     "low": rate.low, "high": rate.high})  # fmt: skip
    return pd.DataFrame(rows).sort_values(["engine", "prep_agree"], ignore_index=True)


def subset_table(per_complex: pd.DataFrame) -> pd.DataFrame:
    """
    How many arms the signal needs: AUROC over every subset of k primary arms, k = 1..5.

    For each subset, agreement is counted over its arms only. Also reported, for each k, is
    the success rate of baseline poses that every arm in the subset agrees with, averaged
    over subsets, and the share of complexes that receive that unanimous verdict.

    Args:
        per_complex (pd.DataFrame): Output of `collect`, with one `d_<arm>` column per arm.

    Returns:
        pd.DataFrame: One row per (engine, k): number of subsets, mean, min and max AUROC,
            the best subset, and the mean unanimous success rate and coverage.
    """
    rows = []
    for e, g in per_complex.groupby("engine", sort=False):
        agree = {a: (g[f"d_{a}"] < SUCCESS_THRESHOLD).astype(int) for a in PRIMARY_ARMS}
        for k in range(1, len(PRIMARY_ARMS) + 1):
            scores, unanimous_rate, unanimous_share = {}, [], []
            for subset in itertools.combinations(PRIMARY_ARMS, k):
                count = sum(agree[a] for a in subset)
                scores[subset] = auroc(count, g["success"])
                unanimous = count == k
                unanimous_rate.append(g.loc[unanimous, "success"].mean())
                unanimous_share.append(unanimous.mean())
            best = max(scores, key=scores.get)
            rows.append({
                "engine": e,
                "k": k,
                "subsets": len(scores),
                "auroc_mean": float(np.mean(list(scores.values()))),
                "auroc_min": min(scores.values()),
                "auroc_max": scores[best],
                "best_subset": "+".join(best),
                "unanimous_success": float(np.nanmean(unanimous_rate)),
                "unanimous_share": float(np.mean(unanimous_share)),
            })  # fmt: skip
    return pd.DataFrame(rows)


def consensus_table(per_complex: pd.DataFrame) -> pd.DataFrame:
    """
    Consensus pose against the baseline pose, per engine and metric, with paired tests.

    The four engines form one Holm family per metric.

    Args:
        per_complex (pd.DataFrame): Output of `collect`.

    Returns:
        pd.DataFrame: One row per (metric, engine): rates, difference, discordant counts,
            exact McNemar p and Holm-adjusted p.
    """
    rows = []
    for metric, (base_col, cons_col) in {
        "top1": ("success", "consensus_success"),
        "pb_valid": ("pb_success", "consensus_pb_success"),
    }.items():
        block = []
        for e, g in per_complex.groupby("engine", sort=False):
            g = g.set_index("id")
            c = compare(g[base_col], g[cons_col])
            block.append({
                "metric": metric,
                "engine": e,
                "n": len(g),
                "baseline": g[base_col].mean(),
                "consensus": g[cons_col].mean(),
                "diff": g[cons_col].mean() - g[base_col].mean(),
                "gained": len(c.only_b),
                "lost": len(c.only_a),
                "p": c.p_value,
                "changed_pose": float((g["consensus_arm"] != BASELINE).mean()),
            })  # fmt: skip
        for row, p in zip(block, holm_adjust([r["p"] for r in block]), strict=True):
            row["p_holm"] = p
        rows += block
    return pd.DataFrame(rows)


def oracle_table(results: str | Path) -> pd.DataFrame:
    """
    Upper bound for choosing the preparation per complex: success if any candidate succeeds.

    Compared with the same bound over seed replicates, where they exist, to separate the
    headroom from preparation from that of rerunning.

    Args:
        results (str | Path): Root containing `<engine>/<arm>/` run directories.

    Returns:
        pd.DataFrame: One row per engine: baseline, best-of-arms and best-of-seeds success.
    """
    rows = []
    for e, r in discover(results).items():
        base = successes(r.baseline)
        arms = pd.concat(
            [base, *(successes(f).reindex(base.index, fill_value=False) for f in r.arms.values())],
            axis=1,
        )
        row = {
            "engine": e,
            "n": len(base),
            "baseline": base.mean(),
            "best_of_arms": arms.any(axis=1).mean(),
            "n_arms": len(r.arms),
        }
        if r.seeds:
            seeds = pd.concat(
                [
                    base,
                    *(successes(f).reindex(base.index, fill_value=False) for f in r.seeds.values()),
                ],
                axis=1,
            )
            row |= {"best_of_seeds": seeds.any(axis=1).mean(), "n_seeds": len(r.seeds)}
        rows.append(row)
    return pd.DataFrame(rows)


def main(argv: Sequence[str] | None = None) -> None:
    """
    Run the full agreement analysis, write its tables as CSV and, optionally, its figure.

    Args:
        argv (Sequence[str] | None): Command-line arguments; sys.argv if None.
    """
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("results", type=Path, help="root with <engine>/<arm>/ runs")
    parser.add_argument("out", type=Path, help="directory for the CSV tables")
    parser.add_argument("--figures", type=Path, help="directory for the agreement figure")
    args = parser.parse_args(argv)
    RDLogger.DisableLog("rdApp.*")
    args.out.mkdir(parents=True, exist_ok=True)
    per_complex = collect(args.results)
    tables = {
        "agreement_per_complex": per_complex,
        "agreement_signals": signal_table(per_complex),
        "agreement_levels": level_table(per_complex),
        "agreement_subsets": subset_table(per_complex),
        "agreement_consensus": consensus_table(per_complex),
        "agreement_oracle": oracle_table(args.results),
    }
    for name, table in tables.items():
        table.to_csv(args.out / f"{name}.csv", index=False)
        print(
            f"\n== {name}\n{table.round(3).to_string(index=False)}"
            if name != "agreement_per_complex"
            else ""
        )
    if args.figures:
        # Imported here: plotting needs matplotlib, which the tables don't.
        from dockprepbench.analysis.figures import plot_agreement_levels

        for path in plot_agreement_levels(
            tables["agreement_levels"], args.figures / "pose_agreement"
        ):
            print(f"Wrote {path}")


if __name__ == "__main__":
    main()
