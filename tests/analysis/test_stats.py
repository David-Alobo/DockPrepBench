import math

import pandas as pd
import pytest

from dockprepbench.analysis.stats import (
    compare,
    holm_adjust,
    load_run,
    success_rate,
    successes,
    wilson_interval,
)


def flags(**hits):
    return pd.Series(hits, dtype=bool)


@pytest.fixture
def run_csv(tmp_path):
    path = tmp_path / "results.csv"
    pd.DataFrame(
        {
            "id": ["A", "B", "C", "D"],
            "top1_rmsd": [0.5, 1.5, 3.0, math.nan],
            "top1_score": [-9.0, -8.0, -7.0, math.nan],
            "top1_pb_valid": [True, False, True, False],
            "top1_failed_checks": ["", "bond_lengths", "", ""],
            "n_poses": [9, 9, 9, 0],
            "seconds": [1.0, 1.0, 1.0, 1.0],
            "error": ["", "", "", "RuntimeError: boom"],
        }
    ).to_csv(path, index=False)
    return path


def test_load_run_keeps_blank_text_as_empty(run_csv):
    run = load_run(run_csv)

    assert run.loc["A", "error"] == ""
    assert run.loc["D", "error"] == "RuntimeError: boom"
    assert run["top1_pb_valid"].dtype == bool


def test_successes_counts_failures_as_misses(run_csv):
    run = load_run(run_csv)

    assert successes(run).to_dict() == {"A": True, "B": True, "C": False, "D": False}
    assert successes(run, pb_valid=True).to_dict() == {
        "A": True,
        "B": False,
        "C": False,
        "D": False,
    }


def test_wilson_interval_matches_reference_values():
    # Reference: Wilson interval for 8/10 at 95% is about [0.490, 0.943].
    low, high = wilson_interval(8, 10)

    assert low == pytest.approx(0.490, abs=1e-3)
    assert high == pytest.approx(0.943, abs=1e-3)


@pytest.mark.parametrize(("k", "n"), [(0, 5), (5, 5)])
def test_wilson_interval_stays_in_unit_range(k, n):
    low, high = wilson_interval(k, n)

    assert 0.0 <= low <= k / n <= high <= 1.0


def test_wilson_interval_rejects_bad_counts():
    with pytest.raises(ValueError):
        wilson_interval(3, 0)


def test_success_rate_reports_counts():
    rate = success_rate(flags(A=True, B=False, C=True, D=True))

    assert (rate.successes, rate.total) == (3, 4)
    assert rate.value == pytest.approx(0.75)
    assert "75.0%" in str(rate)


def test_compare_pairs_on_shared_ids():
    a = flags(A=True, B=True, C=False, D=False, E=True)
    b = flags(A=True, B=False, C=True, D=False, F=True)

    result = compare(a, b)

    assert result.n == 4  # E and F are not shared
    assert result.both == 1
    assert result.only_a == ["B"]
    assert result.only_b == ["C"]
    assert result.neither == 1
    assert result.p_value == pytest.approx(1.0)


def test_compare_detects_one_sided_change():
    ids = [f"X{i}" for i in range(12)]
    a = pd.Series(False, index=ids)
    b = pd.Series(True, index=ids)

    result = compare(a, b)

    # 12 discordant pairs all favouring b: two-sided exact p = 2 * 0.5**12
    assert result.p_value == pytest.approx(2 * 0.5**12)


def test_compare_identical_runs_has_p_one():
    a = flags(A=True, B=False)

    assert compare(a, a).p_value == 1.0


def test_compare_requires_overlap():
    with pytest.raises(ValueError):
        compare(flags(A=True), flags(B=True))


def test_holm_adjust_matches_hand_calculation():
    # Sorted: 0.005*4=0.02, 0.01*3=0.03, 0.03*2=0.06, 0.04*1=0.04 -> monotone max 0.06
    adjusted = holm_adjust([0.01, 0.04, 0.03, 0.005])

    assert adjusted == pytest.approx([0.03, 0.06, 0.06, 0.02])


def test_holm_adjust_caps_at_one_and_keeps_single_p():
    assert holm_adjust([0.6, 0.7]) == pytest.approx([1.0, 1.0])
    assert holm_adjust([0.04]) == pytest.approx([0.04])


def test_paired_difference_interval_contains_estimate_and_is_symmetric_at_zero():
    from dockprepbench.analysis.stats import paired_difference_interval

    low, high = paired_difference_interval(both=180, only_a=18, only_b=48, n=428)
    assert low < 30 / 428 < high
    assert 0.0 < low  # 48 vs 18 discordant is a clear improvement
    low0, high0 = paired_difference_interval(both=100, only_a=0, only_b=0, n=428)
    assert low0 == pytest.approx(-high0)
    with pytest.raises(ValueError):
        paired_difference_interval(both=400, only_a=20, only_b=20, n=428)
