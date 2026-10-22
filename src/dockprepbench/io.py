"""Reading and writing molecular structure files."""

from pathlib import Path

from Bio.PDB import PDBParser
from Bio.PDB.Structure import Structure


def load_pdb_file(file_path: str | Path) -> Structure:
    """
    Load and parse a PDB file using Biopython.

    Args:
        file_path (str | Path): Path to the .pdb file.

    Returns:
        Structure: A Biopython Structure object representing the parsed PDB.

    Raises:
        FileNotFoundError: If the file does not exist.
        ValueError: If the file contains no atoms.
    """
    path = Path(file_path)
    if not path.is_file():
        raise FileNotFoundError(f"PDB file not found: {path}")

    parser = PDBParser(QUIET=True)  # QUIET suppresses warnings
    structure = parser.get_structure(path.stem, str(path))
    if not any(True for _ in structure.get_atoms()):
        raise ValueError(f"No atoms parsed from PDB file: {path}")
    return structure
