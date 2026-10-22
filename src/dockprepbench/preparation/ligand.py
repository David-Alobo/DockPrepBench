"""Ligand preparation: protonation state and starting conformer."""

from dimorphite_dl import protonate_smiles
from rdkit import Chem
from rdkit.Chem import AllChem

from dockprepbench.chem import skeleton

PROTONATION_CHOICES = ("supplied", "dimorphite")
CONFORMER_CHOICES = ("supplied", "rdkit")


def prepare_ligand(
    ligand: Chem.Mol, protonation: str, conformer: str, ph: float, seed: int = 0
) -> Chem.Mol:
    """
    Prepare a ligand for docking according to a preparation protocol.

    Protonation is applied first, so a regenerated conformer is built for the
    protonated molecule.

    Args:
        ligand (Chem.Mol): Ligand with explicit hydrogens and a 3D conformer.
        protonation (str): "supplied" keeps the input state; "dimorphite" assigns the
            dominant state at `ph` with Dimorphite-DL.
        conformer (str): "supplied" keeps the input heavy-atom geometry; "rdkit" embeds
            a new conformer (ETKDGv3, then MMFF) from `seed`.
        ph (float): pH for protonation.
        seed (int): Random seed for conformer generation.

    Returns:
        Chem.Mol: Prepared ligand with explicit hydrogens and one conformer.

    Raises:
        ValueError: If a choice is not recognised.
    """
    if protonation not in PROTONATION_CHOICES:
        raise ValueError(f"protonation must be one of {PROTONATION_CHOICES}, got {protonation!r}")
    if conformer not in CONFORMER_CHOICES:
        raise ValueError(f"conformer must be one of {CONFORMER_CHOICES}, got {conformer!r}")

    if protonation == "dimorphite":
        ligand = protonate_ligand(ligand, ph)
    if conformer == "rdkit":
        ligand = generate_conformer(ligand, seed)
    return ligand


def protonate_ligand(ligand: Chem.Mol, ph: float) -> Chem.Mol:
    """
    Set a ligand's protonation state to the dominant one at a pH, keeping its geometry.

    Dimorphite-DL works on SMILES, so its molecule supplies the chemistry (charges,
    hydrogen counts and bond orders, which can shift in e.g. amidinium groups) and the
    input supplies heavy-atom coordinates, matched atom by atom on the charge-free
    skeleton. Hydrogens are then rebuilt around the heavy atoms.

    Args:
        ligand (Chem.Mol): Ligand with a 3D conformer.
        ph (float): pH at which to assign the protonation state.

    Returns:
        Chem.Mol: Protonated ligand with explicit hydrogens.

    Raises:
        ValueError: If the protonated molecule cannot be mapped back onto the input.
    """
    heavy = Chem.RemoveAllHs(ligand)
    (state,) = protonate_smiles(
        Chem.MolToSmiles(heavy), ph_min=ph, ph_max=ph, precision=0.0, max_variants=1
    )
    protonated = Chem.MolFromSmiles(state)
    match = skeleton(heavy).GetSubstructMatch(skeleton(protonated))
    if len(match) != heavy.GetNumAtoms():
        raise ValueError(f"Cannot map protonated state {state} onto the input ligand")

    positions = heavy.GetConformer().GetPositions()
    conformer = Chem.Conformer(protonated.GetNumAtoms())
    for protonated_idx, heavy_idx in enumerate(match):
        conformer.SetAtomPosition(protonated_idx, positions[heavy_idx].tolist())
    protonated.AddConformer(conformer, assignId=True)
    Chem.AssignStereochemistryFrom3D(protonated)
    return Chem.AddHs(protonated, addCoords=True)


def generate_conformer(ligand: Chem.Mol, seed: int = 0) -> Chem.Mol:
    """
    Replace a ligand's conformer with a new ETKDGv3 embedding, minimised with MMFF94.

    Args:
        ligand (Chem.Mol): Ligand with explicit hydrogens.
        seed (int): Random seed, so the same seed always gives the same conformer.

    Returns:
        Chem.Mol: Copy of the ligand with a single new conformer.

    Raises:
        ValueError: If embedding fails.
    """
    mol = Chem.Mol(ligand)
    mol.RemoveAllConformers()
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    if AllChem.EmbedMolecule(mol, params) != 0:
        raise ValueError(f"Conformer embedding failed for {Chem.MolToSmiles(mol)}")
    if AllChem.MMFFHasAllMoleculeParams(mol):
        AllChem.MMFFOptimizeMolecule(mol)
    else:
        AllChem.UFFOptimizeMolecule(mol)
    return mol
