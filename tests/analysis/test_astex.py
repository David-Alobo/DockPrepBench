import numpy as np
import pandas as pd

from dockprepbench.analysis.agreement import PRIMARY_ARMS
from dockprepbench.analysis.astex import MIN_NEGATIVES, h2_h3_table, h4_table


def per_complex(n=80, seed=2, engine="vina"):
    rng = np.random.default_rng(seed)
    success = rng.random(n) < 0.6
    frame = pd.DataFrame({
        "engine": engine,
        "id": [f"C{i}" for i in range(n)],
        "success": success,
        "score": -7.0 + rng.normal(0, 1, n),
        "cnn": np.nan,
    })  # fmt: skip
    for arm in PRIMARY_ARMS:
        frame[f"d_{arm}"] = np.where(rng.random(n) < np.where(success, 0.9, 0.2), 1.0, 5.0)
    frame["prep_agree"] = sum((frame[f"d_{a}"] < 2).astype(int) for a in PRIMARY_ARMS)
    return frame


def test_h2_h3_supports_a_strong_signal_and_skips_too_few_negatives():
    t = h2_h3_table(per_complex()).set_index("hypothesis")
    assert t.loc["H2", "tested"] and t.loc["H2", "supported"]
    assert t.loc["H3", "auroc"] > t.loc["H3", "score_auroc"]

    easy = per_complex()
    easy["success"] = True
    easy.loc[: MIN_NEGATIVES - 2, "success"] = False  # fewer incorrect poses than the minimum
    t = h2_h3_table(easy)
    assert not t["tested"].any() and not t["supported"].any()


def test_h4_compares_unanimous_success_with_posebusters():
    frame = per_complex()
    levels = pd.DataFrame(
        [{"engine": "vina", "prep_agree": 5, "n": 166, "rate": 0.837, "low": 0.774, "high": 0.886}]
    )
    row = h4_table(frame, levels).iloc[0]
    assert row["n"] == (frame["prep_agree"] == 5).sum()
    assert row["pb_rate"] == 0.837
    assert isinstance(row["consistent"], (bool, np.bool_))
