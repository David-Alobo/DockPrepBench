"""Common interface shared by all docking engines."""

from abc import ABC, abstractmethod
from dataclasses import dataclass
from pathlib import Path
from typing import ClassVar

from rdkit import Chem


@dataclass(frozen=True)
class Box:
    """Axis-aligned search box in Å, in the receptor frame."""

    center: tuple[float, float, float]
    size: tuple[float, float, float]


@dataclass(frozen=True)
class DockingResult:
    """Docked poses ranked best-first, with one score per pose (lower is better)."""

    poses: list[Chem.Mol]
    scores: list[float]

    def __post_init__(self) -> None:
        if len(self.poses) != len(self.scores):
            raise ValueError(f"Got {len(self.poses)} poses but {len(self.scores)} scores")


class DockingEngine(ABC):
    """A docking program that places one ligand into one receptor."""

    name: ClassVar[str]

    @abstractmethod
    def dock(self, receptor: Path, ligand: Chem.Mol, box: Box) -> DockingResult:
        """
        Dock a ligand into a receptor.

        Args:
            receptor (Path): Prepared, protonated receptor PDB file.
            ligand (Chem.Mol): Ligand with explicit hydrogens and a 3D conformer.
            box (Box): Search space.

        Returns:
            DockingResult: Ranked poses and their scores.
        """
