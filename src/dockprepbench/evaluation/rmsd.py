"""Ligand pose RMSD and top-1 docking success."""

import math
from collections.abc import Iterable

from rdkit import Chem
from rdkit.Chem import rdMolAlign

SUCCESS_THRESHOLD = 2.0  # Å, the conventional cutoff for a correct pose


def pose_rmsd(pose: Chem.Mol, reference: Chem.Mol) -> float:
    """
    Compute the symmetry-corrected heavy-atom RMSD between a docked pose and a reference.

    The pose is not superimposed onto the reference: docked poses already share the
    receptor frame, so any offset is part of the error being measured.

    Args:
        pose (Chem.Mol): Docked ligand with 3D coordinates.
        reference (Chem.Mol): Crystal ligand with 3D coordinates.

    Returns:
        float: RMSD in Å, minimised over symmetry-equivalent atom mappings.

    Raises:
        RuntimeError: If the two molecules do not share the same heavy-atom graph.
    """
    return rdMolAlign.CalcRMS(Chem.RemoveAllHs(pose), Chem.RemoveAllHs(reference))


def is_success(rmsd: float, threshold: float = SUCCESS_THRESHOLD) -> bool:
    """Return True if a pose RMSD is below the success threshold. NaN counts as failure."""
    return not math.isnan(rmsd) and rmsd < threshold


def top1_success_rate(rmsds: Iterable[float], threshold: float = SUCCESS_THRESHOLD) -> float:
    """
    Compute the fraction of complexes whose top-ranked pose is a success.

    Args:
        rmsds (Iterable[float]): Top-1 pose RMSD per complex. Use NaN for failed runs.
        threshold (float): Success cutoff in Å.

    Returns:
        float: Success rate between 0 and 1.

    Raises:
        ValueError: If no RMSDs are given.
    """
    outcomes = [is_success(r, threshold) for r in rmsds]
    if not outcomes:
        raise ValueError("Cannot compute a success rate from no RMSDs")
    return sum(outcomes) / len(outcomes)
