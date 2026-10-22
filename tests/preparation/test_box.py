import pytest
from rdkit import Chem
from rdkit.Geometry import Point3D

from dockprepbench.preparation.box import box_from_ligand


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
