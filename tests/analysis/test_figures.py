import math

import pandas as pd
import pytest

from dockprepbench.analysis.figures import (
    SERIES,
    _spread,
    plot_rmsd_cdf,
    plot_success_rates,
    summary_table,
)


def run(rmsds, valid=None):
    n = len(rmsds)
    return pd.DataFrame(
        {
            "top1_rmsd": rmsds,
            "top1_pb_valid": valid if valid is not None else [True] * n,
            "error": ["" if not math.isnan(r) else "RuntimeError: x" for r in rmsds],
        },
        index=[f"C{i}" for i in range(n)],
    )


@pytest.fixture
def runs():
    return {
        "baseline": run([0.5, 1.0, 3.0, math.nan], [True, False, True, False]),
        "box_padding_2": run([0.5, 1.5, 1.8, 4.0]),
    }


def test_summary_table_rates_and_failures(runs):
    table = summary_table(runs)

    assert list(table.index) == ["baseline", "box_padding_2"]
    assert table.loc["baseline", "top1"] == pytest.approx(0.5)
    assert table.loc["baseline", "pb_valid"] == pytest.approx(0.25)
    assert table.loc["baseline", "failed"] == 1
    assert table.loc["box_padding_2", "top1"] == pytest.approx(0.75)
    assert (table["top1_low"] <= table["top1"]).all() and (
        table["top1"] <= table["top1_high"]
    ).all()


def test_plots_write_every_format(runs, tmp_path):
    written = plot_success_rates(summary_table(runs), tmp_path / "rates")
    written += plot_rmsd_cdf(runs, tmp_path / "cdf")

    assert sorted(p.name for p in written) == ["cdf.pdf", "cdf.png", "rates.pdf", "rates.png"]
    assert all(p.stat().st_size > 0 for p in written)


def test_rmsd_cdf_refuses_more_runs_than_colours(tmp_path):
    many = {f"r{i}": run([1.0]) for i in range(len(SERIES) + 1)}

    with pytest.raises(ValueError, match="At most"):
        plot_rmsd_cdf(many, tmp_path / "cdf")


def test_spread_keeps_order_and_minimum_gap():
    values = [80.0, 81.0, 50.0, 80.5]

    placed = _spread(values, min_gap=5.0)

    assert [anchor for anchor, _ in placed] == values
    positions = sorted(label for _, label in placed)
    assert all(b - a >= 5.0 - 1e-9 for a, b in zip(positions, positions[1:], strict=False))
    assert placed[2][1] < placed[0][1]  # the low curve's label stays below the others


def test_spread_leaves_separated_labels_alone():
    assert _spread([10.0, 50.0], min_gap=5.0) == [(10.0, 10.0), (50.0, 50.0)]
