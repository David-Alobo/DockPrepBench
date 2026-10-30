import pandas as pd
from test_report import results  # noqa: F401 - fixture, same directory; tests are not a package

from dockprepbench.analysis.latex import _p, _signed, arm_table, sensitivity_table, write_latex
from dockprepbench.analysis.report import write_report


def test_p_values_keep_two_significant_figures():
    assert _p(0.2) == "0.20"
    assert _p(0.0234) == "0.023"
    assert _p(0.999) == "1.0"
    assert _p(4.1e-4) == r"$4.1\times10^{-4}$"
    assert _p(float("nan")) == "--"


def test_changes_carry_a_sign():
    assert _signed(6.78) == "+6.8"
    assert _signed(-1.9) == "$-1.9$"


def test_arm_table_lists_baseline_then_arms_with_their_comparisons(results, tmp_path):  # noqa: F811
    write_report(results, tmp_path)
    body = arm_table(
        pd.read_csv(tmp_path / "summary.csv"), pd.read_csv(tmp_path / "comparisons.csv"), "top1"
    )
    rows = [line for line in body.splitlines() if line.endswith(r"\\")][1:]  # skip the header

    # Engines follow the figures' order: Vina before Gnina.
    assert rows[0].startswith("Vina & Baseline")
    assert "& -- & -- & -- & --" in rows[0]
    box2 = next(r for r in rows if "Box padding 2" in r)
    assert "+50.0 [+5.3, +73.3] & 3 & 0" in box2  # 33.3% -> 83.3%; 3 gained, none lost


def test_write_latex_wraps_whole_environments(results, tmp_path):  # noqa: F811
    write_report(results, tmp_path)

    top1 = (tmp_path / "top1_table.tex").read_text(encoding="utf-8")
    assert top1.startswith(r"\begin{tabular}{llrrrrr}") and top1.rstrip().endswith(r"\end{tabular}")


def test_write_latex_skips_tables_whose_csv_is_missing(tmp_path):
    assert write_latex(tmp_path) == []


def test_sensitivity_table_is_a_longtable_body():
    frame = pd.DataFrame(
        [
            {"engine": "vina", "arm": "box_padding_2", "property": "heavy_atoms",
             "flipped_value": 29, "unflipped_value": 24, "effect": 0.21,
             "effect_kind": "rank-biserial r", "p_holm": 0.057},
            {"engine": "vina", "arm": "box_padding_2", "property": "cofactor",
             "flipped_value": 0.5, "unflipped_value": 0.6, "effect": 0.76,
             "effect_kind": "odds ratio", "p_holm": 1.0},
        ]
    )  # fmt: skip

    body = sensitivity_table(frame)

    assert r"\endhead" in body and r"\endfoot" in body
    assert r"Heavy atoms & 29 & 24 & $r = 0.21$ & 0.057" in body
    assert r"50.0\% & 60.0\% & OR 0.76 & 1.0" in body
