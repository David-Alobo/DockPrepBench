"""Conversion between RDKit/PDB and the PDBQT format used by AutoDock-family engines."""

from meeko import (
    MoleculePreparation,
    PDBQTMolecule,
    PDBQTWriterLegacy,
    Polymer,
    RDKitMolCreate,
    ResidueChemTemplates,
)
from rdkit import Chem


def receptor_to_pdbqt(pdb_text: str) -> str:
    """
    Convert a protonated receptor PDB to rigid-receptor PDBQT.

    Args:
        pdb_text (str): PDB contents. Every residue must match a meeko template,
            so hydrogens must already be present.

    Returns:
        str: Rigid receptor PDBQT.

    Raises:
        meeko.PolymerCreationError: If a residue cannot be matched to a template.
    """
    polymer = Polymer.from_pdb_string(
        pdb_text, ResidueChemTemplates.create_from_defaults(), MoleculePreparation()
    )
    rigid, _flexible = PDBQTWriterLegacy.write_from_polymer(polymer)
    return rigid


def ligand_to_pdbqt(ligand: Chem.Mol) -> str:
    """
    Convert a ligand to PDBQT, keeping its SMILES so bond orders survive docking.

    Args:
        ligand (Chem.Mol): Ligand with explicit hydrogens and a 3D conformer.

    Returns:
        str: Ligand PDBQT with torsion tree.

    Raises:
        ValueError: If meeko cannot write the ligand.
    """
    setups = MoleculePreparation().prepare(ligand)
    if len(setups) != 1:
        raise ValueError(f"Expected one molecule setup, got {len(setups)}")
    pdbqt, ok, error = PDBQTWriterLegacy.write_string(setups[0])
    if not ok:
        raise ValueError(f"meeko failed to write ligand PDBQT: {error}")
    return pdbqt


def pdbqt_to_poses(pdbqt: str) -> list[Chem.Mol]:
    """
    Convert docked PDBQT output back to RDKit molecules.

    Args:
        pdbqt (str): Docking output written from a meeko-prepared ligand,
            one MODEL per pose.

    Returns:
        list[Chem.Mol]: One molecule per pose, in file order, each with a single conformer.
    """
    (mol,) = RDKitMolCreate.from_pdbqt_mol(PDBQTMolecule(pdbqt, skip_typing=True))
    return [Chem.Mol(mol, confId=conf.GetId()) for conf in mol.GetConformers()]
