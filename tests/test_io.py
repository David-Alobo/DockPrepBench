import pytest

from dockprepbench.io import load_pdb_file


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
