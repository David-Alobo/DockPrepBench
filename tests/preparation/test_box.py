import pytest
from rdkit import Chem
from rdkit.Geometry import Point3D

from dockprepbench.preparation.box import (
    GYRATION_FACTOR,
    box_from_gyration,
    box_from_ligand,
    radius_of_gyration,
)


def ligand_at(*points):
    mol = Chem.RWMol(Chem.MolFromSmiles("C" * len(points)))
    conf = Chem.Conformer(len(points))
    for i, p in enumerate(points):
        conf.SetAtomPosition(i, Point3D(*p))
    mol.AddConformer(conf)
    return mol.GetMol()


def test_box_is_centred_on_extent_with_padding():
    mol = ligand_at((0, 0, 0), (2, 4, 6), (1, 1, 1))

    box = box_from_ligand(mol, padding=4.0)

    assert box.center == pytest.approx((1.0, 2.0, 3.0))
    assert box.size == pytest.approx((10.0, 12.0, 14.0))


def test_box_ignores_hydrogens():
    mol = Chem.AddHs(ligand_at((0, 0, 0), (2, 0, 0)), addCoords=True)

    box = box_from_ligand(mol, padding=0.0)

    assert box.size == pytest.approx((2.0, 0.0, 0.0))


def test_radius_of_gyration_of_two_atoms_is_half_their_distance():
    assert radius_of_gyration(ligand_at((0, 0, 0), (4, 0, 0))) == pytest.approx(2.0)


def test_radius_of_gyration_ignores_hydrogens():
    mol = Chem.AddHs(ligand_at((0, 0, 0), (4, 0, 0)), addCoords=True)

    assert radius_of_gyration(mol) == pytest.approx(2.0)


def test_gyration_box_is_a_cube_scaled_to_the_docked_ligand():
    reference = ligand_at((0, 0, 0), (2, 0, 0), (0, 6, 0))
    docked = ligand_at((10, 10, 10), (14, 10, 10))  # radius of gyration 2 Å, anywhere

    box = box_from_gyration(reference, docked)

    edge = 2.0 * GYRATION_FACTOR
    assert box.size == pytest.approx((edge, edge, edge))
    assert box.center == pytest.approx((2 / 3, 2.0, 0.0))  # the reference's centroid
