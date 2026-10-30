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


def holm_adjust(p_values: list[float]) -> list[float]:
    """
    Adjust p-values for multiple comparisons with Holm's step-down method.

    Controls the family-wise error rate (the chance of any false positive across all
    comparisons) without assuming independence, and is uniformly more powerful than
    Bonferroni.

    Args:
        p_values (list[float]): Raw p-values, one per comparison.

    Returns:
        list[float]: Adjusted p-values, in the input order.
    """
    m = len(p_values)
    order = sorted(range(m), key=p_values.__getitem__)
    adjusted = [0.0] * m
    running_max = 0.0
    for rank, i in enumerate(order):
        running_max = max(running_max, min(1.0, (m - rank) * p_values[i]))
        adjusted[i] = running_max
    return adjusted


def paired_difference_interval(
    both: int, only_a: int, only_b: int, n: int, z: float = Z_95
) -> tuple[float, float]:
    """
    Confidence interval for the difference of two paired proportions, rate(b) - rate(a).

    Newcombe's hybrid score interval (Newcombe 1998, method 10), recommended by Fagerland
    et al. (2014): it combines the Wilson intervals of the two rates with their observed
    correlation, so unlike two separate Wilson intervals it reflects the pairing.

    Args:
        both (int): Complexes that succeed in both runs.
        only_a (int): Complexes that succeed in run a only.
        only_b (int): Complexes that succeed in run b only.
        n (int): Complexes compared.
        z (float): Normal quantile for the confidence level.

    Returns:
        tuple[float, float]: Lower and upper bounds of rate(b) - rate(a).

    Raises:
        ValueError: If the counts are inconsistent with n.
    """
    neither = n - both - only_a - only_b
    if n <= 0 or min(both, only_a, only_b, neither) < 0:
        raise ValueError(f"Inconsistent counts: {both=}, {only_a=}, {only_b=}, {n=}")
    p_a, p_b = (both + only_a) / n, (both + only_b) / n
    l_a, u_a = wilson_interval(both + only_a, n, z)
    l_b, u_b = wilson_interval(both + only_b, n, z)
    margins = (both + only_a) * (only_b + neither) * (both + only_b) * (only_a + neither)
    phi = (both * neither - only_a * only_b) / math.sqrt(margins) if margins else 0.0
    d = p_b - p_a
    low = d - math.sqrt(
        max((p_b - l_b) ** 2 - 2 * phi * (p_b - l_b) * (u_a - p_a) + (u_a - p_a) ** 2, 0.0)
    )
    high = d + math.sqrt(
        max((u_b - p_b) ** 2 - 2 * phi * (u_b - p_b) * (p_a - l_a) + (p_a - l_a) ** 2, 0.0)
    )
    return max(-1.0, low), min(1.0, high)


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
