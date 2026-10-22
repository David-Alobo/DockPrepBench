import pytest

PDB_TEXT = """\
ATOM      1  N   ALA A   1      11.104   6.134  -6.504  1.00  0.00           N
ATOM      2  CA  ALA A   1      11.639   6.071  -5.147  1.00  0.00           C
ATOM      3  C   ALA A   1      13.149   5.879  -5.168  1.00  0.00           C
ATOM      4  O   ALA A   1      13.741   5.394  -6.135  1.00  0.00           O
END
"""


@pytest.fixture
def pdb_file(tmp_path):
    """A minimal four-atom PDB file named 1abc.pdb."""
    path = tmp_path / "1abc.pdb"
    path.write_text(PDB_TEXT)
    return path


@pytest.fixture
def peptide_pdb(tmp_path):
    """A protonated, zwitterionic AGSA tetrapeptide that meeko can template-match."""
    from rdkit import Chem
    from rdkit.Chem import AllChem

    pep = Chem.MolFromSequence("AGSA")
    for atom in pep.GetAtoms():
        info = atom.GetPDBResidueInfo()
        if info.GetName().strip() == "OXT":
            atom.SetFormalCharge(-1)
        elif info.GetName().strip() == "N" and info.GetResidueNumber() == 1:
            atom.SetFormalCharge(1)
    Chem.SanitizeMol(pep)
    pep = Chem.AddHs(pep, addResidueInfo=True)
    AllChem.EmbedMolecule(pep, randomSeed=1)
    # meeko needs each residue's atoms to be contiguous; AddHs appends hydrogens at the end.
    order = sorted(
        range(pep.GetNumAtoms()),
        key=lambda i: (pep.GetAtomWithIdx(i).GetPDBResidueInfo().GetResidueNumber(), i),
    )
    path = tmp_path / "peptide.pdb"
    path.write_text(Chem.MolToPDBBlock(Chem.RenumberAtoms(pep, order)))
    return path
