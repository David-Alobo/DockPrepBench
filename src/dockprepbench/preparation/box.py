"""Search box definition."""

from rdkit import Chem

from dockprepbench.engines.base import Box

DEFAULT_PADDING = 4.0  # Å added on each side, matching smina's --autobox_add default


def box_from_ligand(ligand: Chem.Mol, padding: float = DEFAULT_PADDING) -> Box:
    """
    Define a search box around a ligand's heavy atoms.

    Args:
        ligand (Chem.Mol): Ligand with a 3D conformer, typically the crystal ligand.
        padding (float): Margin in Å added on each side of the ligand's extent.

    Returns:
        Box: Box centred on the ligand's extent.
    """
    coords = Chem.RemoveAllHs(ligand).GetConformer().GetPositions()
    low, high = coords.min(axis=0), coords.max(axis=0)
    center = (low + high) / 2
    size = high - low + 2 * padding
    return Box(center=tuple(center.tolist()), size=tuple(size.tolist()))
