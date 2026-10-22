import math

import pytest
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Geometry import Point3D

from dockprepbench.evaluation.rmsd import is_success, pose_rmsd, top1_success_rate


def embedded(smiles: str) -> Chem.Mol:
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    AllChem.EmbedMolecule(mol, randomSeed=42)
    return mol


def translated(mol: Chem.Mol, dx: float) -> Chem.Mol:
    moved = Chem.Mol(mol)
    conf = moved.GetConformer()
    for i in range(moved.GetNumAtoms()):
        p = conf.GetAtomPosition(i)
        conf.SetAtomPosition(i, Point3D(p.x + dx, p.y, p.z))
    return moved


def test_identical_pose_has_zero_rmsd():
    mol = embedded("c1ccccc1C(=O)N")

    assert pose_rmsd(mol, mol) == pytest.approx(0.0)


def test_translation_is_not_superimposed_away():
    mol = embedded("c1ccccc1C(=O)N")

    assert pose_rmsd(translated(mol, 3.0), mol) == pytest.approx(3.0)


def test_symmetric_atom_swap_has_zero_rmsd():
    # Isobutane: swapping two equivalent methyl carbons changes indices, not geometry.
    ref = Chem.RemoveHs(embedded("CC(C)C"))
    swapped = Chem.Mol(ref)
    conf = swapped.GetConformer()
    p0, p2 = conf.GetAtomPosition(0), conf.GetAtomPosition(2)
    conf.SetAtomPosition(0, p2)
    conf.SetAtomPosition(2, p0)

    assert pose_rmsd(swapped, ref) == pytest.approx(0.0, abs=1e-6)


def test_hydrogens_are_ignored():
    mol = embedded("CCO")

    assert pose_rmsd(mol, Chem.RemoveHs(mol)) == pytest.approx(0.0)


def test_different_molecules_raise():
    with pytest.raises(RuntimeError):
        pose_rmsd(embedded("CCO"), embedded("CCN"))


@pytest.mark.parametrize(
    ("rmsd", "expected"), [(0.5, True), (1.99, True), (2.0, False), (5.0, False), (math.nan, False)]
)
def test_is_success(rmsd, expected):
    assert is_success(rmsd) is expected


def test_top1_success_rate_counts_nan_as_failure():
    assert top1_success_rate([1.0, 3.0, math.nan, 0.5]) == pytest.approx(0.5)


def test_top1_success_rate_custom_threshold():
    assert top1_success_rate([1.0, 2.5], threshold=3.0) == pytest.approx(1.0)


def test_top1_success_rate_empty_raises():
    with pytest.raises(ValueError):
        top1_success_rate([])
