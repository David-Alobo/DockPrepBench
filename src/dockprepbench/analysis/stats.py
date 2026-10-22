"""Success rates with confidence intervals, and paired comparison of two runs."""

import math
from dataclasses import dataclass
from pathlib import Path

import pandas as pd
from scipy.stats import binomtest

from dockprepbench.evaluation.rmsd import SUCCESS_THRESHOLD, is_success

Z_95 = 1.959963984540054  # two-sided 95% normal quantile


@dataclass(frozen=True)
class Rate:
    """A success proportion with its 95% Wilson score interval."""

    successes: int
    total: int
    low: float
    high: float

    @property
    def value(self) -> float:
        """Point estimate of the proportion."""
        return self.successes / self.total

    def __str__(self) -> str:
        return (
            f"{100 * self.value:.1f}% [{100 * self.low:.1f}, {100 * self.high:.1f}] "
            f"({self.successes}/{self.total})"
        )


@dataclass(frozen=True)
class Comparison:
    """Paired outcome of two runs on the same complexes."""

    n: int
    both: int
    only_a: list[str]  # ids that succeed in run a but not in run b
    only_b: list[str]  # ids that succeed in run b but not in run a
    p_value: float  # exact McNemar test on the discordant pairs

    @property
    def neither(self) -> int:
        """Number of complexes that fail in both runs."""
        return self.n - self.both - len(self.only_a) - len(self.only_b)


def load_run(path: str | Path) -> pd.DataFrame:
    """
    Read a benchmark results CSV, keeping blank errors and failed checks as empty strings.

    Args:
        path (str | Path): CSV written by `dockprepbench run`.

    Returns:
        pd.DataFrame: One row per complex, indexed by complex id.
    """
    df = pd.read_csv(path)
    for column in ("top1_failed_checks", "error"):
        df[column] = df[column].fillna("")
    return df.set_index("id")


def successes(
    run: pd.DataFrame, pb_valid: bool = False, threshold: float = SUCCESS_THRESHOLD
) -> pd.Series:
    """
    Flag each complex as a top-1 success. Failed runs (NaN RMSD) count as misses.

    Args:
        run (pd.DataFrame): Output of `load_run`.
        pb_valid (bool): Also require the top pose to pass every PoseBusters check.
        threshold (float): RMSD cutoff in Å.

    Returns:
        pd.Series: Boolean per complex id.
    """
    hit = run["top1_rmsd"].map(lambda r: is_success(r, threshold))
    if pb_valid:
        hit &= run["top1_pb_valid"]
    return hit


def wilson_interval(k: int, n: int, z: float = Z_95) -> tuple[float, float]:
    """
    Compute the Wilson score interval for a binomial proportion.

    Unlike the normal approximation it stays inside [0, 1] and behaves well near 0 and 1.

    Args:
        k (int): Number of successes.
        n (int): Number of trials.
        z (float): Normal quantile for the confidence level.

    Returns:
        tuple[float, float]: Lower and upper bounds.

    Raises:
        ValueError: If n is not positive or k is outside [0, n].
    """
    if n <= 0 or not 0 <= k <= n:
        raise ValueError(f"Need 0 <= k <= n and n > 0, got k={k}, n={n}")
    p = k / n
    centre = (p + z**2 / (2 * n)) / (1 + z**2 / n)
    half = z * math.sqrt(p * (1 - p) / n + z**2 / (4 * n**2)) / (1 + z**2 / n)
    return max(0.0, centre - half), min(1.0, centre + half)


def success_rate(hits: pd.Series) -> Rate:
    """
    Summarise success flags as a proportion with a 95% confidence interval.

    Args:
        hits (pd.Series): Boolean success flag per complex.

    Returns:
        Rate: Proportion and Wilson interval.
    """
    k, n = int(hits.sum()), len(hits)
    low, high = wilson_interval(k, n)
    return Rate(successes=k, total=n, low=low, high=high)


def compare(hits_a: pd.Series, hits_b: pd.Series) -> Comparison:
    """
    Compare two runs on the complexes they share, with an exact McNemar test.

    Only discordant complexes (success in exactly one run) carry information about a
    difference; under no difference each is equally likely to favour either run.

    Args:
        hits_a (pd.Series): Success flags of run a, indexed by complex id.
        hits_b (pd.Series): Success flags of run b, indexed by complex id.

    Returns:
        Comparison: Counts, the ids that changed outcome, and the two-sided p-value.

    Raises:
        ValueError: If the runs share no complexes.
    """
    a, b = hits_a.align(hits_b, join="inner")
    if a.empty:
        raise ValueError("The runs share no complexes")
    only_a = sorted(a.index[a & ~b])
    only_b = sorted(a.index[~a & b])
    discordant = len(only_a) + len(only_b)
    p_value = binomtest(len(only_b), discordant, 0.5).pvalue if discordant else 1.0
    return Comparison(
        n=len(a), both=int((a & b).sum()), only_a=only_a, only_b=only_b, p_value=p_value
    )
