import math

import pandas as pd
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem
from test_report import IDS, write_run  # same directory; tests are not a package

from dockprepbench.analysis.report import discover
from dockprepbench.analysis.sensitivity import (
    cofactors,
    ligand_properties,
    seed_sensitive,
    sensitivity_table,
)
from dockprepbench.datasets.posebusters import Complex


def make_complex(tmp_path, cid, smiles, hetatms=()):
    d = tmp_path / cid
    d.mkdir()
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    AllChem.EmbedMolecule(mol, randomSeed=1)
    for name in ("ligand", "start"):
        with Chem.SDWriter(str(d / f"{name}.sdf")) as w:
            w.write(mol)
    lines = ["ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C"]
    lines += [
        f"HETATM    2  C1  {h:>3} A 101       1.000   1.000   1.000  1.00  0.00           C"
        for h in hetatms
    ]
    (d / "protein.pdb").write_text("\n".join(lines) + "\nEND\n")
    return Complex(
        id=cid,
        protein=d / "protein.pdb",
        crystal_ligand=d / "ligand.sdf",
        start_ligand=d / "start.sdf",
    )


def test_cofactors_skip_additives_and_the_ligand_itself(tmp_path):
    entry = make_complex(tmp_path, "1ABC_LIG", "CCO", hetatms=("HEM", "HOH", "GOL", "LIG", "NAD"))

    assert cofactors(entry) == ["HEM", "NAD"]


def test_ligand_properties(tmp_path):
    acid = make_complex(tmp_path, "1AAA_ACE", "CC(=O)O", hetatms=("HEM",))  # acetate at pH 7.4
    benzene = make_complex(tmp_path, "2BBB_BNZ", "c1ccccc1")

    props = ligand_properties([acid, benzene])

    assert props.loc["1AAA_ACE", "charge_changed"]
    assert not props.loc["2BBB_BNZ", "charge_changed"]
    assert props.loc["1AAA_ACE", "cofactor"] and not props.loc["2BBB_BNZ", "cofactor"]
    assert props.loc["2BBB_BNZ", ["heavy_atoms", "rotatable_bonds", "abs_charge"]].tolist() == [
        6,
        0,
        0,
    ]


@pytest.fixture
def runs(tmp_path):
    #                         C0   C1   C2   C3   C4   C5
    write_run(tmp_path, "vina", "baseline", [1.0, 1.0, 1.0, 3.0, 3.0, 3.0])
    write_run(tmp_path, "vina", "baseline_seed1", [1.0, 3.0, 1.0, 3.0, 3.0, 3.0])  # C1 noisy
    write_run(tmp_path, "vina", "box_padding_2", [1.0, 3.0, 3.0, 1.0, 3.0, 3.0])  # C1-C3 flip
    write_run(tmp_path, "smina", "baseline", [1.0, 1.0, 1.0, 3.0, 3.0, 3.0])
    write_run(tmp_path, "smina", "box_padding_2", [1.0] * 6)  # C3-C5 flip
    return discover(tmp_path)


@pytest.fixture
def props():
    return pd.DataFrame(
        {
            "heavy_atoms": [10, 30, 32, 34, 12, 14],
            "rotatable_bonds": [1, 2, 3, 4, 5, 6],
            "abs_charge": [0, 0, 1, 0, 0, 1],
            "charge_changed": [True, True, False, None, False, False],
            "cofactor": [False, True, True, True, False, False],
        },
        index=IDS,
    )


def test_seed_sensitive_marks_complexes_that_differ_between_replicates(runs):
    assert seed_sensitive(runs["vina"]).tolist() == [False, True, False, False, False, False]
    assert seed_sensitive(runs["smina"]) is None


def test_sensitivity_compares_flipped_with_unflipped(runs, props):
    t = sensitivity_table(runs, props).set_index(["engine", "arm", "property"])
    vina = t.loc[("vina", "box_padding_2")]

    assert set(vina["n_flipped"]) == {3}
    # Flipped C1-C3 have 30, 32, 34 heavy atoms; the others 10, 12, 14.
    assert vina.loc["heavy_atoms", ["flipped_value", "unflipped_value"]].tolist() == [32, 12]
    assert vina.loc["heavy_atoms", "effect"] == pytest.approx(1.0)  # complete separation
    # Cofactor: all 3 flipped have one, none of the 3 others.
    assert vina.loc["cofactor", ["flip_rate_with", "flip_rate_without"]].tolist() == [1.0, 0.0]
    assert math.isinf(vina.loc["cofactor", "effect"])


def test_unknown_values_are_left_out(runs, props):
    t = sensitivity_table(runs, props).set_index(["engine", "arm", "property"])

    # charge_changed is unknown for C3, so 2 of the 3 flipped complexes count.
    row = t.loc[("vina", "box_padding_2", "charge_changed")]
    assert row["flip_rate_with"] == pytest.approx(1 / 2)  # C0 unflipped, C1 flipped


def test_seed_sensitivity_is_tested_only_with_replicates_and_p_is_holm_adjusted(runs, props):
    t = sensitivity_table(runs, props)

    assert "seed_sensitive" in set(t[t.engine == "vina"].property)
    assert "seed_sensitive" not in set(t[t.engine == "smina"].property)
    assert (t["p_holm"] >= t["p"] - 1e-12).all()
