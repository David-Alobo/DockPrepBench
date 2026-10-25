"""Search box definition."""

import numpy as np
from rdkit import Chem

from dockprepbench.engines.base import Box

DEFAULT_PADDING = 4.0  # Å added on each side, matching smina's --autobox_add default
# How the box is sized: a padding around the reference ligand's extent, or a cube scaled to
# the docked ligand's radius of gyration.
BOX_SIZE_CHOICES = ("padding", "gyration")
# Edge of the cubic box per Å of radius of gyration, the optimum Feinstein and Brylinski
# (J. Cheminform. 2015, 7, 18) found for AutoDock Vina.
GYRATION_FACTOR = 2.857


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


def radius_of_gyration(ligand: Chem.Mol) -> float:
    """
    Radius of gyration of a ligand's heavy atoms, unweighted by mass.

    Args:
        ligand (Chem.Mol): Ligand with a 3D conformer.

    Returns:
        float: Root-mean-square distance of the heavy atoms from their centroid, in Å.
    """
    coords = Chem.RemoveAllHs(ligand).GetConformer().GetPositions()
    return float(np.sqrt(((coords - coords.mean(axis=0)) ** 2).sum(axis=1).mean()))


def box_from_gyration(reference: Chem.Mol, ligand: Chem.Mol) -> Box:
    """
    Define a cubic search box scaled to the docked ligand's size.

    Following Feinstein and Brylinski (2015), the edge is 2.857 times the radius of gyration
    of the ligand's heavy atoms, computed on the conformer that is docked (here a low-energy
    start conformer, as in their protocol), and the box is centred on the geometric centre
    of the reference ligand.

    Args:
        reference (Chem.Mol): Ligand that marks the site, typically the crystal ligand.
        ligand (Chem.Mol): Ligand to be docked, with a 3D conformer.

    Returns:
        Box: Cubic box at the reference ligand's centroid.
    """
    heavy = Chem.RemoveAllHs(reference)
    center = heavy.GetConformer().GetPositions().mean(axis=0)
    edge = GYRATION_FACTOR * radius_of_gyration(ligand)
    return Box(center=tuple(center.tolist()), size=(edge, edge, edge))
