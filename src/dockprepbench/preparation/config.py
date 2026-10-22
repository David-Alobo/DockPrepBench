"""Named preparation protocols: the variables this benchmark varies."""

import tomllib
from dataclasses import asdict, dataclass, fields
from pathlib import Path
from typing import Any

from dockprepbench.preparation.box import DEFAULT_PADDING
from dockprepbench.preparation.ligand import CONFORMER_CHOICES, PROTONATION_CHOICES
from dockprepbench.preparation.receptor import DEFAULT_PH


@dataclass(frozen=True)
class PrepConfig:
    """
    One preparation protocol. Every field is a choice whose effect on docking is measured.

    Attributes:
        name (str): Short identifier, used in result paths.
        ph (float): pH for receptor protonation.
        box_padding (float): Å added on each side of the crystal ligand's extent.
        ligand_protonation (str): "supplied" or "dimorphite" (dominant state at `ph`).
        ligand_conformer (str): "supplied" start conformer, or "rdkit" to regenerate it.
    """

    name: str = "baseline"
    ph: float = DEFAULT_PH
    box_padding: float = DEFAULT_PADDING
    ligand_protonation: str = "supplied"
    ligand_conformer: str = "supplied"

    def __post_init__(self) -> None:
        if self.ligand_protonation not in PROTONATION_CHOICES:
            raise ValueError(
                f"ligand_protonation must be one of {PROTONATION_CHOICES}, "
                f"got {self.ligand_protonation!r}"
            )
        if self.ligand_conformer not in CONFORMER_CHOICES:
            raise ValueError(
                f"ligand_conformer must be one of {CONFORMER_CHOICES}, "
                f"got {self.ligand_conformer!r}"
            )

    def receptor_key(self) -> str:
        """
        Identify the settings that change the prepared receptor, for caching.

        Configs with the same key share prepared receptors, so varying only the box
        does not repeat the slow receptor preparation.

        Returns:
            str: Cache key such as "ph7.4".
        """
        return f"ph{self.ph:g}"

    def to_dict(self) -> dict[str, Any]:
        """
        Return the config as a plain dict, for writing next to results.

        Returns:
            dict[str, Any]: Field names and values.
        """
        return asdict(self)


def load_config(path: str | Path) -> PrepConfig:
    """
    Read a preparation config from a TOML file.

    Fields left out of the file keep their baseline defaults, so a variant file only
    needs to state what it changes. The name defaults to the file's stem.

    Args:
        path (str | Path): TOML file with top-level keys matching PrepConfig fields.

    Returns:
        PrepConfig: The parsed config.

    Raises:
        ValueError: If the file contains a key that is not a PrepConfig field.
    """
    path = Path(path)
    values = tomllib.loads(path.read_text(encoding="utf-8"))
    unknown = set(values) - {f.name for f in fields(PrepConfig)}
    if unknown:
        raise ValueError(f"Unknown config keys in {path}: {sorted(unknown)}")
    values.setdefault("name", path.stem)
    return PrepConfig(**values)
