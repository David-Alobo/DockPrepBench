"""Reading and writing molecular structure files."""

from pathlib import Path

from Bio.PDB import PDBParser
from Bio.PDB.Structure import Structure
from rdkit import Chem


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


def load_sdf_file(file_path: str | Path) -> list[Chem.Mol]:
    """
    Load all molecules from an SDF file, keeping explicit hydrogens.

    Args:
        file_path (str | Path): Path to the .sdf file.

    Returns:
        list[Chem.Mol]: Molecules in file order (e.g. docked poses ranked by score).

    Raises:
        FileNotFoundError: If the file does not exist.
        ValueError: If any record fails to parse, or the file contains no molecules.
    """
    path = Path(file_path)
    if not path.is_file():
        raise FileNotFoundError(f"SDF file not found: {path}")

    try:
        mols = list(Chem.SDMolSupplier(str(path), removeHs=False))
    except OSError as err:  # RDKit rejects empty files as invalid input
        raise ValueError(f"No molecules in SDF file: {path}") from err
    if not mols:
        raise ValueError(f"No molecules in SDF file: {path}")
    if any(mol is None for mol in mols):
        raise ValueError(f"Failed to parse a record in SDF file: {path}")
    return mols


def write_sdf_file(mols: list[Chem.Mol], file_path: str | Path) -> Path:
    """
    Write molecules to an SDF file atomically, keeping their properties and hydrogens.

    The file is written under a temporary name and renamed, so an interrupted write
    never leaves a partial file behind.

    Args:
        mols (list[Chem.Mol]): Molecules to write, in order.
        file_path (str | Path): Destination .sdf file; parent directories are created.

    Returns:
        Path: The written file.
    """
    path = Path(file_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".sdf.tmp")
    with Chem.SDWriter(str(tmp)) as writer:
        for mol in mols:
            writer.write(mol)
    tmp.replace(path)
    return path
