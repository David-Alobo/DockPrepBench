import pandas as pd
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from dockprepbench.analysis.followup import PADDINGS, _inside, interaction_table
from dockprepbench.engines.base import Box


def write_run(root, name, outcomes):
    folder = root / "vina" / name
    folder.mkdir(parents=True)
    rows = [{"id": cid, "top1_rmsd": 1.0 if ok else 5.0, "top1_score": -7.0,
             "top1_pb_valid": True, "top1_failed_checks": "", "top1_centroid_shift": 0.5,
             "best_rmsd": 1.0, "best_rank": 1, "n_poses": 9, "seconds": 1.0, "error": ""}
            for cid, ok in outcomes.items()]  # fmt: skip
    pd.DataFrame(rows).to_csv(folder / "results.csv", index=False)


def test_interaction_compares_protonation_at_each_padding(tmp_path):
    ids = [f"C{i}" for i in range(20)]
    for sup, prot in PADDINGS.values():
        write_run(tmp_path, sup, {c: i < 10 for i, c in enumerate(ids)})
        write_run(tmp_path, prot, {c: i < 14 for i, c in enumerate(ids)})  # 4 gained
    table = interaction_table(tmp_path / "vina").set_index("padding")
    assert list(table.index) == [2, 4, 6, 8]
    assert (table["gained"] == 4).all() and (table["lost"] == 0).all()
    assert table.loc[2, "diff_pts"] == pytest.approx(20.0)
    assert table.loc[4, "vs_baseline_padding_p"] != table.loc[4, "vs_baseline_padding_p"]  # NaN
    assert (table["p_holm"] >= table["p"]).all()


def test_inside_checks_every_heavy_atom():
    mol = Chem.AddHs(Chem.MolFromSmiles("CCO"))
    AllChem.EmbedMolecule(mol, randomSeed=1)
    centre = tuple(Chem.RemoveAllHs(mol).GetConformer().GetPositions().mean(axis=0))
    assert _inside(mol, Box(center=centre, size=(20.0, 20.0, 20.0)))
    assert not _inside(mol, Box(center=centre, size=(0.5, 0.5, 0.5)))
