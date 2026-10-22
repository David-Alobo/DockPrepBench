"""Physical validity of docked poses with the PoseBusters checks."""

import warnings
from pathlib import Path

from posebusters import PoseBusters
from rdkit import Chem

# PoseBusters' own RMSD check; dropped because RMSD is scored separately (evaluation.rmsd).
RMSD_CHECK = "rmsd_≤_2å"


def validity_checks(pose: Chem.Mol, crystal: Chem.Mol, protein: str | Path) -> dict[str, bool]:
    """
    Run the PoseBusters redocking checks on a pose, except RMSD.

    These cover chemical consistency with the crystal ligand (formula, bonds,
    stereochemistry), internal geometry (bond lengths and angles, ring planarity, clashes,
    energy) and contacts with the protein, its cofactors and waters.

    Args:
        pose (Chem.Mol): Docked ligand pose.
        crystal (Chem.Mol): Crystal ligand, used for the chemistry checks.
        protein (str | Path): Protein PDB to check contacts against. Use the original
            structure, with cofactors, to match the published PoseBusters protocol.

    Returns:
        dict[str, bool]: Pass/fail per check, keyed by PoseBusters check name.
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")  # PoseBusters warns per check on routine inputs
        report = PoseBusters(config="redock").bust([pose], crystal, str(protein))
    checks = report.iloc[0].drop(RMSD_CHECK)
    return {name: bool(passed) for name, passed in checks.items()}


def is_pb_valid(checks: dict[str, bool]) -> bool:
    """
    Return True if a pose passes every PoseBusters check.

    Args:
        checks (dict[str, bool]): Output of `validity_checks`.

    Returns:
        bool: Whether the pose is PB-valid.
    """
    return all(checks.values())
