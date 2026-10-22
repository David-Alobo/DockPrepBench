import pytest

from dockprepbench.engines.pdbqt import receptor_to_pdbqt
from dockprepbench.preparation.receptor import protonate_receptor

WATER = "HETATM  999  O   HOH W 101      30.000  30.000  30.000  1.00  0.00           O\n"


def is_hydrogen(line):
    return line.startswith("ATOM") and line[76:78].strip() == "H"


def heavy_atoms_only(pdb_text):
    lines = pdb_text.splitlines(keepends=True)
    return "".join(line for line in lines if not is_hydrogen(line))


def test_protonate_receptor_adds_hydrogens_and_drops_heterogens(peptide_pdb, tmp_path):
    bare = tmp_path / "bare.pdb"
    bare.write_text(heavy_atoms_only(peptide_pdb.read_text()).replace("END", WATER + "END"))
    assert not any(is_hydrogen(line) for line in bare.read_text().splitlines())

    prepared = protonate_receptor(bare)

    assert any(is_hydrogen(line) for line in prepared.splitlines())
    assert "HOH" not in prepared
    assert receptor_to_pdbqt(prepared)  # meeko accepts the result


def atom_coords(pdb_text):
    return {
        (line[21], line[22:27].strip(), line[12:16].strip()): line[30:54]
        for line in pdb_text.splitlines()
        if line.startswith("ATOM")
    }


def test_protonate_receptor_keeps_crystal_coordinates(peptide_pdb, tmp_path):
    # Drop the serine OG so PDBFixer has a heavy atom to rebuild and relax.
    text = heavy_atoms_only(peptide_pdb.read_text())
    text = "".join(line for line in text.splitlines(True) if line[12:20] != " OG  SER")
    bare = tmp_path / "bare.pdb"
    bare.write_text(text)

    prepared = atom_coords(protonate_receptor(bare))

    original = atom_coords(text)
    assert ("A", "3", "OG") in prepared
    for key, coords in original.items():
        assert [float(c) for c in prepared[key].split()] == pytest.approx(
            [float(c) for c in coords.split()], abs=1e-3
        )


def test_protonate_receptor_is_reproducible(peptide_pdb, tmp_path):
    # Rebuilding a side chain exercises every step that involves random placement.
    text = heavy_atoms_only(peptide_pdb.read_text())
    text = "".join(line for line in text.splitlines(True) if line[12:20] != " OG  SER")
    bare = tmp_path / "bare.pdb"
    bare.write_text(text)

    assert protonate_receptor(bare) == protonate_receptor(bare)


def test_protonate_receptor_removes_isolated_residues(peptide_pdb, tmp_path):
    lone_ala = (
        "ATOM    900  N   ALA X   1      20.000  20.000  20.000  1.00  0.00           N\n"
        "ATOM    901  CA  ALA X   1      21.458  20.000  20.000  1.00  0.00           C\n"
        "ATOM    902  C   ALA X   1      22.009  21.420  20.000  1.00  0.00           C\n"
        "ATOM    903  O   ALA X   1      21.236  22.381  20.000  1.00  0.00           O\n"
        "ATOM    904  CB  ALA X   1      21.986  19.215  21.194  1.00  0.00           C\n"
        "TER\n"
    )
    bare = tmp_path / "bare.pdb"
    bare.write_text(heavy_atoms_only(peptide_pdb.read_text()).replace("END", lone_ala + "END"))

    prepared = protonate_receptor(bare)

    assert not any(line[21] == "X" for line in prepared.splitlines() if line.startswith("ATOM"))
    assert any(line[21] == "A" for line in prepared.splitlines() if line.startswith("ATOM"))
