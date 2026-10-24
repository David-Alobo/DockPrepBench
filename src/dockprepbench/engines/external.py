"""Helpers for engines that run as external programs."""

import shutil
import subprocess
from collections.abc import Sequence
from pathlib import Path

from rdkit import Chem

DEFAULT_TIMEOUT = 3600  # seconds per complex; guards against a hung process


class EngineError(RuntimeError):
    """An external docking program failed or produced no usable output."""


def resolve_command(command: Sequence[str] | None, program: str) -> list[str]:
    """
    Return the command prefix for an external program.

    Args:
        command (Sequence[str] | None): Explicit prefix, e.g. a path to the binary or
            `[python, fake_engine.py]` in tests. If None, `program` is looked up on PATH.
        program (str): Executable name to look up.

    Returns:
        list[str]: Command prefix.

    Raises:
        EngineError: If the program is not on PATH.
    """
    if command:
        return list(command)
    path = shutil.which(program)
    if path is None:
        raise EngineError(f"{program} not found on PATH; install it or pass its path explicitly")
    return [path]


def run(
    args: list[str], cwd: Path, timeout: int = DEFAULT_TIMEOUT, env: dict | None = None
) -> None:
    """
    Run an external program, raising with its output on failure.

    Args:
        args (list[str]): Full command line.
        cwd (Path): Working directory.
        timeout (int): Seconds before the process is killed.
        env (dict | None): Environment for the process; None inherits the current one.

    Raises:
        EngineError: If the program exits non-zero or times out.
    """
    try:
        proc = subprocess.run(
            args, cwd=cwd, capture_output=True, text=True, timeout=timeout, env=env, check=False
        )
    except subprocess.TimeoutExpired as err:
        raise EngineError(f"{Path(args[0]).name} timed out after {timeout} s") from err
    if proc.returncode != 0:
        tail = (proc.stderr or proc.stdout).strip().splitlines()[-5:]
        raise EngineError(f"{Path(args[0]).name} exited {proc.returncode}: {' | '.join(tail)}")


def write_ligand(ligand: Chem.Mol, path: Path) -> Path:
    """
    Write a ligand (with hydrogens and coordinates) to an SDF file.

    Args:
        ligand (Chem.Mol): Prepared ligand.
        path (Path): Destination .sdf file.

    Returns:
        Path: The written file.
    """
    with Chem.SDWriter(str(path)) as writer:
        writer.write(ligand)
    return path


def read_poses(path: Path) -> list[Chem.Mol]:
    """
    Read docked poses from an SDF file in file order, keeping hydrogens and properties.

    Args:
        path (Path): Output .sdf from the engine.

    Returns:
        list[Chem.Mol]: Poses that RDKit could parse.

    Raises:
        EngineError: If the file is missing or holds no readable pose.
    """
    if not path.is_file():
        raise EngineError(f"No output file {path.name}")
    poses = [m for m in Chem.SDMolSupplier(str(path), removeHs=False) if m is not None]
    if not poses:
        raise EngineError(f"No readable poses in {path.name}")
    return poses
