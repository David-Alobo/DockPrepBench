import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from dockprepbench.engines.pdbqt import ligand_to_pdbqt, pdbqt_to_poses, receptor_to_pdbqt
from dockprepbench.evaluation.rmsd import pose_rmsd


@pytest.fixture
def benzamide():
    mol = Chem.AddHs(Chem.MolFromSmiles("NC(=O)c1ccccc1"))
    AllChem.EmbedMolecule(mol, randomSeed=1)
    return mol


def test_ligand_round_trip_keeps_chemistry_and_coordinates(benzamide):
    (pose,) = pdbqt_to_poses(ligand_to_pdbqt(benzamide))

    assert Chem.MolToSmiles(Chem.RemoveHs(pose)) == Chem.MolToSmiles(Chem.RemoveHs(benzamide))
    assert pose_rmsd(pose, benzamide) == pytest.approx(0.0, abs=1e-2)


def test_multi_model_pdbqt_gives_one_pose_per_model(benzamide):
    pdbqt = ligand_to_pdbqt(benzamide)
    models = "".join(f"MODEL {i}\n{pdbqt}ENDMDL\n" for i in (1, 2, 3))

    poses = pdbqt_to_poses(models)

    assert len(poses) == 3
    assert all(p.GetNumConformers() == 1 for p in poses)


def test_ligand_without_hydrogens_is_rejected():
    mol = Chem.MolFromSmiles("CCO")
    AllChem.Compute2DCoords(mol)

    with pytest.raises(ValueError, match="explicit Hs"):
        ligand_to_pdbqt(mol)


def test_receptor_to_pdbqt_types_atoms(peptide_pdb):
    pdbqt = receptor_to_pdbqt(peptide_pdb.read_text())

    atom_lines = [line for line in pdbqt.splitlines() if line.startswith("ATOM")]
    assert atom_lines
    assert {line[77:79].strip() for line in atom_lines} >= {"C", "OA", "N", "HD"}
