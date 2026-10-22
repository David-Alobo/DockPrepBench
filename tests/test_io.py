import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from dockprepbench.io import load_pdb_file, load_sdf_file


def test_load_pdb_file_uses_stem_as_id(pdb_file):
    structure = load_pdb_file(pdb_file)

    assert structure.id == "1abc"
    assert len(list(structure.get_atoms())) == 4


def test_load_pdb_file_missing_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        load_pdb_file(tmp_path / "missing.pdb")


def test_load_pdb_file_empty_raises(tmp_path):
    pdb = tmp_path / "empty.pdb"
    pdb.write_text("END\n")

    with pytest.raises(ValueError):
        load_pdb_file(pdb)


def test_load_sdf_file_keeps_order_and_hydrogens(tmp_path):
    sdf = tmp_path / "poses.sdf"
    writer = Chem.SDWriter(str(sdf))
    for smiles in ("CCO", "CCN"):
        mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
        AllChem.EmbedMolecule(mol, randomSeed=0)
        writer.write(mol)
    writer.close()

    mols = load_sdf_file(sdf)

    assert [Chem.MolToSmiles(Chem.RemoveHs(m)) for m in mols] == ["CCO", "CCN"]
    assert mols[0].GetNumAtoms() == 9


def test_load_sdf_file_missing_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        load_sdf_file(tmp_path / "missing.sdf")


def test_load_sdf_file_empty_raises(tmp_path):
    sdf = tmp_path / "empty.sdf"
    sdf.write_text("")

    with pytest.raises(ValueError):
        load_sdf_file(sdf)
