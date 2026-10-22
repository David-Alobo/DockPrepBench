import numpy as np
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolTransforms
from rdkit.Geometry import Point3D

from dockprepbench.evaluation.validity import RMSD_CHECK, is_pb_valid, validity_checks
from dockprepbench.io import load_pdb_file


@pytest.fixture
def peptide_center(peptide_pdb):
    return np.mean([a.coord for a in load_pdb_file(peptide_pdb).get_atoms()], axis=0)


def placed(mol, center):
    """Copy of mol translated so its heavy-atom centroid sits at center."""
    moved = Chem.Mol(mol)
    conf = moved.GetConformer()
    shift = center - Chem.RemoveHs(mol).GetConformer().GetPositions().mean(axis=0)
    for i in range(moved.GetNumAtoms()):
        conf.SetAtomPosition(i, Point3D(*(np.array(conf.GetAtomPosition(i)) + shift)))
    return moved


@pytest.fixture
def ethanol():
    mol = Chem.AddHs(Chem.MolFromSmiles("CCO"))
    AllChem.EmbedMolecule(mol, randomSeed=1)
    AllChem.MMFFOptimizeMolecule(mol)
    return mol


def test_rmsd_check_is_excluded(ethanol, peptide_pdb, peptide_center):
    pose = placed(ethanol, peptide_center + [0, 0, 6])

    checks = validity_checks(pose, ethanol, peptide_pdb)

    assert RMSD_CHECK not in checks
    assert "bond_lengths" in checks and "minimum_distance_to_protein" in checks


def test_clash_with_protein_is_flagged(ethanol, peptide_pdb, peptide_center):
    pose = placed(ethanol, peptide_center)  # sits on top of the peptide atoms

    checks = validity_checks(pose, ethanol, peptide_pdb)

    assert not checks["minimum_distance_to_protein"]
    assert not is_pb_valid(checks)


def test_stretched_bond_is_flagged(ethanol, peptide_pdb, peptide_center):
    pose = placed(ethanol, peptide_center + [0, 0, 6])
    oxygen = next(a for a in pose.GetAtoms() if a.GetSymbol() == "O")
    carbon = oxygen.GetNeighbors()[0]
    rdMolTransforms.SetBondLength(pose.GetConformer(), carbon.GetIdx(), oxygen.GetIdx(), 2.5)

    checks = validity_checks(pose, ethanol, peptide_pdb)

    assert not checks["bond_lengths"]


def test_is_pb_valid():
    assert is_pb_valid({"a": True, "b": True})
    assert not is_pb_valid({"a": True, "b": False})
