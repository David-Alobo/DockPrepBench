import math

import pandas as pd
import pytest

from dockprepbench.analysis.report import (
    comparison_table,
    discover,
    flips_table,
    noise_table,
    outcome_class,
    read_ids,
    restrict,
    summary_table,
    transition_table,
    write_report,
)

IDS = [f"C{i}" for i in range(6)]


def write_run(root, engine, arm, rmsd, valid=None, best=None, shift=None, error=None):
    n = len(rmsd)
    d = root / engine / arm
    d.mkdir(parents=True)
    pd.DataFrame(
        {
            "id": IDS[:n],
            "top1_rmsd": rmsd,
            "top1_score": -8.0,
            "top1_pb_valid": valid or [True] * n,
            "top1_failed_checks": "",
            "top1_centroid_shift": shift or [1.0] * n,
            "best_rmsd": best or rmsd,
            "best_rank": 1,
            "n_poses": 9,
            "seconds": 1.0,
            "error": error or [""] * n,
        }
    ).to_csv(d / "results.csv", index=False)


@pytest.fixture
def results(tmp_path):
    #                 C0   C1   C2   C3   C4   C5
    write_run(
        tmp_path,
        "vina",
        "baseline",
        [1.0, 1.0, 3.0, 3.0, 3.0, 3.0],
        best=[1.0, 1.0, 1.5, 3.0, 3.0, 3.0],
    )
    write_run(tmp_path, "vina", "baseline_seed1", [1.0, 3.0, 3.0, 3.0, 3.0, 3.0])
    write_run(tmp_path, "vina", "baseline_seed2", [1.0, 1.0, 1.0, 3.0, 3.0, 3.0])
    write_run(tmp_path, "vina", "box_padding_2", [1.0, 1.0, 1.0, 1.0, 1.0, 3.0])
    write_run(tmp_path, "vina", "box_padding_8", [1.0, 3.0, 3.0, 3.0, 3.0, 3.0])
    write_run(tmp_path, "gnina", "baseline", [1.0, 1.0, 1.0, 1.0, 3.0, 3.0])
    write_run(tmp_path, "smina", "box_padding_2", [1.0] * 6)  # no baseline: skipped
    return tmp_path


def test_discover_groups_seeds_and_arms(results):
    runs = discover(results)

    assert list(runs) == ["gnina", "vina"]
    assert sorted(runs["vina"].seeds) == ["baseline_seed1", "baseline_seed2"]
    assert sorted(runs["vina"].arms) == ["box_padding_2", "box_padding_8"]
    assert runs["gnina"].seeds == {} and runs["gnina"].arms == {}


@pytest.mark.parametrize(
    ("row", "expected"),
    [
        ({"error": "RuntimeError: x"}, "failed"),
        ({"top1_rmsd": 1.0, "top1_pb_valid": True}, "valid"),
        ({"top1_rmsd": 1.0, "top1_pb_valid": False}, "invalid"),
        ({"top1_rmsd": 3.0, "best_rmsd": 1.5}, "scoring"),
        ({"top1_rmsd": 3.0, "best_rmsd": 3.0, "top1_centroid_shift": 5.0}, "site"),
        ({"top1_rmsd": 3.0, "best_rmsd": 3.0, "top1_centroid_shift": 2.0}, "mode"),
    ],
)
def test_outcome_class(row, expected):
    full = {"error": "", "top1_rmsd": 3.0, "top1_pb_valid": True, "best_rmsd": 3.0,
            "top1_centroid_shift": 1.0} | row  # fmt: skip
    assert outcome_class(pd.Series(full)) == expected


def test_summary_covers_baselines_arms_and_seeds(results):
    s = summary_table(discover(results)).set_index(["engine", "run"])

    assert len(s) == 6
    assert s.loc[("vina", "baseline"), "top1"] == pytest.approx(2 / 6)
    assert s.loc[("gnina", "baseline"), "top1"] == pytest.approx(4 / 6)


def test_comparisons_are_paired_and_holm_adjusted_per_engine_and_metric(results):
    c = comparison_table(discover(results))
    top1 = c[(c.engine == "vina") & (c.metric == "top1")].set_index("arm")

    assert top1.loc["box_padding_2", ["gained", "lost"]].tolist() == [3, 0]
    assert top1.loc["box_padding_2", "delta_pts"] == pytest.approx(50.0)
    assert (top1["p_holm"] >= top1["p"]).all()
    assert set(top1["family_size"]) == {2}


def test_noise_table_counts_seed_flips(results):
    noise = noise_table(discover(results)).set_index("engine")

    # baseline vs seed1: 1 flip; baseline vs seed2: 1; seed1 vs seed2: 2
    assert noise.loc["vina", "flips_mean"] == pytest.approx(4 / 3)
    assert noise.loc["vina", ["always_succeed", "always_fail", "seed_sensitive"]].tolist() == [
        1,
        3,
        2,
    ]
    assert "gnina" not in noise.index  # no seed replicates


def test_flips_are_compared_with_seed_noise(results):
    runs = discover(results)
    flips = flips_table(runs, noise_table(runs)).set_index("arm")

    assert flips.loc["box_padding_2", "flipped"] == 3
    assert flips.loc["box_padding_2", "vs_noise"] == pytest.approx(3 / (4 / 3))


def test_transitions_count_class_moves(results):
    t = transition_table(discover(results))
    moves = t[(t.arm == "box_padding_2")].set_index(["from", "to"])["n"]

    assert moves[("scoring", "valid")] == 1  # C2: a good pose existed but was not ranked first
    assert moves[("mode", "valid")] == 2  # C3, C4
    assert moves.sum() == 6


def test_write_report_writes_every_table(results, tmp_path):
    written = write_report(results, tmp_path / "tables")

    assert sorted(p.name for p in written) == [
        "comparisons.csv",
        "flips.csv",
        "flips_table.tex",
        "noise.csv",
        "pb_valid_table.tex",
        "summary.csv",
        "top1_table.tex",
        "transitions.csv",
    ]  # fmt: skip (no dataset given, so no sensitivity tables)
    assert not math.isnan(pd.read_csv(tmp_path / "tables" / "noise.csv")["flips_mean"][0])


def test_write_report_needs_a_baseline(tmp_path):
    write_run(tmp_path, "vina", "box_padding_2", [1.0])

    with pytest.raises(ValueError, match="No finished baseline"):
        write_report(tmp_path, tmp_path / "out")


def test_restrict_keeps_only_the_subset_in_every_run(results):
    runs = restrict(discover(results), ["C0", "C2", "C9"])  # C9 is in no run

    assert list(runs["vina"].baseline.index) == ["C0", "C2"]
    assert list(runs["vina"].seeds["baseline_seed1"].index) == ["C0", "C2"]
    assert list(runs["vina"].arms["box_padding_2"].index) == ["C0", "C2"]


def test_write_report_on_a_subset(results, tmp_path):
    write_report(results, tmp_path / "sub", subset=["C2", "C3"])

    summary = pd.read_csv(tmp_path / "sub" / "summary.csv").set_index(["engine", "run"])
    assert summary.loc[("vina", "baseline"), "n"] == 2
    assert summary.loc[("vina", "box_padding_2"), "top1"] == 1.0  # C2 and C3 succeed


def test_read_ids_skips_blanks_and_comments(tmp_path):
    path = tmp_path / "ids.txt"
    path.write_text("# subset\nC0\n\nC2  # note\n")

    assert read_ids(path) == ["C0", "C2"]
