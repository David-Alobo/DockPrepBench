import math

import numpy as np
import pandas as pd
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Geometry import Point3D

from dockprepbench.analysis.agreement import (
    agreement,
    auroc,
    bootstrap_aurocs,
    consensus_index,
    pair_rmsd,
    top_pose,
)


def pose(shift: float = 0.0) -> Chem.Mol:
    mol = Chem.AddHs(Chem.MolFromSmiles("CCOc1ccccc1"))
    AllChem.EmbedMolecule(mol, randomSeed=7)
    conf = mol.GetConformer()
    for i in range(mol.GetNumAtoms()):
        p = conf.GetAtomPosition(i)
        conf.SetAtomPosition(i, Point3D(p.x + shift, p.y, p.z))
    return mol


def test_auroc_perfect_random_and_ties():
    assert auroc([3, 2, 1, 0], [True, True, False, False]) == 1.0
    assert auroc([0, 1, 2, 3], [True, True, False, False]) == 0.0
    assert auroc([1, 1, 1, 1], [True, False, True, False]) == 0.5


def test_auroc_drops_nan_and_needs_both_classes():
    assert auroc([2, float("nan"), 1], [True, True, False]) == 1.0
    assert math.isnan(auroc([1, 2], [True, True]))


def test_agreement_counts_missing_as_disagreement():
    d = pd.DataFrame({"a": [0.5, 3.0], "b": [np.nan, 1.0]})
    assert agreement(d).tolist() == [1, 1]


def test_pair_rmsd_is_shift_without_superposition():
    assert pair_rmsd(pose(), pose(1.5)) == pytest.approx(1.5, abs=1e-6)
    assert math.isnan(pair_rmsd(pose(), None))


def test_consensus_prefers_majority_and_keeps_baseline_on_tie():
    assert consensus_index([pose(5.0), pose(0.0), pose(0.3), pose(0.6)]) == 1
    assert consensus_index([pose(0.0), pose(5.0)]) == 0
    assert consensus_index([None, pose(0.0)]) == 1


def test_bootstrap_differences_are_paired():
    rng = np.random.default_rng(1)
    y = pd.Series(rng.random(200) < 0.5)
    signals = pd.DataFrame({"good": y + rng.normal(0, 0.3, 200), "noise": rng.random(200)})
    t = bootstrap_aurocs(signals, y, reference="noise", n=200)
    assert t.loc["noise", "diff"] == 0
    assert t.loc["good", "low"] <= t.loc["good", "auroc"] <= t.loc["good", "high"]
    assert t.loc["good", "diff_low"] > 0


def test_top_pose_reads_first_pose(tmp_path):
    (tmp_path / "poses").mkdir()
    with Chem.SDWriter(str(tmp_path / "poses" / "X.sdf")) as w:
        w.write(pose(0.0))
        w.write(pose(4.0))
    first = top_pose(tmp_path, "X")
    assert pair_rmsd(first, pose(0.0)) == pytest.approx(0.0, abs=1e-3)
    assert top_pose(tmp_path, "missing") is None


def test_subset_table_counts_subsets_and_full_set_matches():
    from dockprepbench.analysis.agreement import PRIMARY_ARMS, subset_table

    rng = np.random.default_rng(3)
    n = 120
    success = rng.random(n) < 0.5
    frame = pd.DataFrame({"engine": "vina", "success": success})
    for arm in PRIMARY_ARMS:
        # Correct poses mostly agree (small distance), wrong ones mostly don't.
        frame[f"d_{arm}"] = np.where(rng.random(n) < np.where(success, 0.8, 0.3), 1.0, 5.0)
    t = subset_table(frame).set_index("k")
    assert t["subsets"].tolist() == [5, 10, 10, 5, 1]
    full = sum((frame[f"d_{a}"] < 2).astype(int) for a in PRIMARY_ARMS)
    assert t.loc[5, "auroc_mean"] == pytest.approx(auroc(full, frame["success"]))
    assert t.loc[1, "auroc_min"] <= t.loc[1, "auroc_mean"] <= t.loc[1, "auroc_max"]


def test_logistic_fit_recovers_direction_and_validation_runs():
    from dockprepbench.analysis.agreement import (
        PRIMARY_ARMS,
        fit_logistic,
        predict_logistic,
        stratified_folds,
        validation_table,
    )

    rng = np.random.default_rng(5)
    x = rng.normal(size=(300, 1))
    y = (rng.random(300) < 1 / (1 + np.exp(-2 * x[:, 0]))).astype(float)
    beta = fit_logistic(x, y)
    assert beta[1] > 1
    assert auroc(predict_logistic(beta, x), y.astype(bool)) > 0.7

    folds = stratified_folds(y.astype(bool), 5, rng)
    assert set(folds) == set(range(5))
    assert max(abs((y[folds == f]).mean() - y.mean()) for f in range(5)) < 0.05

    success = rng.random(150) < 0.5
    frame = pd.DataFrame(
        {
            "engine": "vina",
            "id": range(150),
            "success": success,
            "score": np.where(success, -8.0, -7.0) + rng.normal(0, 2, 150),
        }
    )
    for arm in PRIMARY_ARMS:
        frame[f"d_{arm}"] = np.where(rng.random(150) < np.where(success, 0.85, 0.3), 1.0, 5.0)
    frame["prep_agree"] = sum((frame[f"d_{a}"] < 2).astype(int) for a in PRIMARY_ARMS)
    table, calibration = validation_table(frame, repeats=2)
    assert set(table["model"]) == {"score", "prep_agree", "score+prep_agree", "best_pair"}
    assert calibration["n"].sum() == 150
