import pandas as pd

from dockprepbench.analysis.report import EngineRuns
from dockprepbench.analysis.robustness import cofactor_failure_table, threshold_table


def run(rows):
    df = pd.DataFrame(
        rows,
        columns=[
            "id",
            "top1_rmsd",
            "best_rmsd",
            "top1_pb_valid",
            "top1_centroid_shift",
            "top1_failed_checks",
            "error",
        ],
    )
    return df.set_index("id")


def test_threshold_moves_misses_between_site_and_mode():
    base = run([("A", 6.0, 5.0, True, 3.0, "", ""), ("B", 7.0, 6.0, True, 9.0, "", "")])
    counts = threshold_table({"vina": EngineRuns("vina", base, {}, {})}, thresholds=(2.0, 4.0))
    by_t = counts.set_index("threshold")
    assert by_t.loc[2.0, "site"] == 2 and by_t.loc[2.0, "mode"] == 0
    assert by_t.loc[4.0, "site"] == 1 and by_t.loc[4.0, "mode"] == 1


def test_cofactor_only_failures_are_counted_for_correct_poses():
    base = run(
        [
            ("A", 1.0, 1.0, False, 0.5, "minimum_distance_to_organic_cofactors", ""),
            ("B", 1.0, 1.0, False, 0.5, "minimum_distance_to_protein", ""),
            ("C", 5.0, 5.0, False, 6.0, "minimum_distance_to_organic_cofactors", ""),
        ]
    )
    table = cofactor_failure_table(
        {"vina": EngineRuns("vina", base, {}, {})}, pd.Series({"A": True, "B": True, "C": True})
    )
    row = table.iloc[0]
    assert row["with_cofactor"] == 3
    assert row["top1_success_with_cofactor"] == 2
    assert row["failed_on_cofactor_only"] == 1
