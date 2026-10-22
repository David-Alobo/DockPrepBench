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

from dockprepbench.engines.base import Box

# Residues that fail template matching are dropped only if they lie this far outside the
# box: Vina's interaction cutoff, so they could not have touched a pose inside it.
BAD_RESIDUE_CLEARANCE = 8.0  # Å


def receptor_to_pdbqt(pdb_text: str, box: Box | None = None) -> str:
    """
    Convert a protonated receptor PDB to rigid-receptor PDBQT.

    Args:
        pdb_text (str): PDB contents. Hydrogens must already be present.
        box (Box | None): Search box. If given, residues that fail template matching are
            dropped (with a logged warning) when they lie more than
            BAD_RESIDUE_CLEARANCE outside it. Without a box, every residue must match.

    Returns:
        str: Rigid receptor PDBQT.

    Raises:
        meeko.PolymerCreationError: If a residue near the box, or any residue when no
            box is given, cannot be matched to a template.
    """
    box_options = {}
    if box is not None:
        box_options = {
            "box_center": list(box.center),
            "box_size": list(box.size),
            "delete_bad_res_from_box_radius": BAD_RESIDUE_CLEARANCE,
        }
    polymer = Polymer.from_pdb_string(
        pdb_text, ResidueChemTemplates.create_from_defaults(), MoleculePreparation(), **box_options
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
