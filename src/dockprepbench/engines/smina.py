"""Smina and Gnina, run as external programs (they share Smina's command line)."""

import os
from collections.abc import Sequence
from pathlib import Path
from tempfile import TemporaryDirectory

from rdkit import Chem

from dockprepbench.engines.base import Box, DockingEngine, DockingResult
from dockprepbench.engines.external import read_poses, resolve_command, run, write_ligand


class SminaEngine(DockingEngine):
    """
    Smina with its default (Vina-derived) scoring function.

    Args:
        exhaustiveness (int): Search effort per run.
        n_poses (int): Maximum number of poses to return.
        seed (int): Random seed. Fixed by default so runs are reproducible.
        cpu (int): Threads per docking run.
        command (Sequence[str] | None): Program to run; default looks up `smina` on PATH.
    """

    name = "smina"
    program = "smina"

    def __init__(
        self,
        exhaustiveness: int = 8,
        n_poses: int = 9,
        seed: int = 42,
        cpu: int = 1,
        command: Sequence[str] | None = None,
    ):
        self.exhaustiveness = exhaustiveness
        self.n_poses = n_poses
        self.seed = seed
        self.cpu = cpu
        self.command = list(command) if command else None

    def extra_args(self) -> list[str]:
        """Engine-specific arguments appended to the shared Smina command line."""
        return []

    def environment(self) -> dict | None:
        """Environment for the docking process; None inherits the current one."""
        return None

    def score_of(self, pose: Chem.Mol) -> float:
        """Score reported for a pose (kcal/mol, lower is better)."""
        return float(pose.GetProp("minimizedAffinity"))

    def dock(self, receptor: Path, ligand: Chem.Mol, box: Box) -> DockingResult:
        with TemporaryDirectory() as tmp:
            work = Path(tmp)
            ligand_file = write_ligand(ligand, work / "ligand.sdf")
            out = work / "poses.sdf"
            (cx, cy, cz), (sx, sy, sz) = box.center, box.size
            args = [
                *resolve_command(self.command, self.program),
                "--receptor", str(Path(receptor).resolve()),
                "--ligand", str(ligand_file),
                "--center_x", f"{cx:.3f}", "--center_y", f"{cy:.3f}", "--center_z", f"{cz:.3f}",
                "--size_x", f"{sx:.3f}", "--size_y", f"{sy:.3f}", "--size_z", f"{sz:.3f}",
                "--exhaustiveness", str(self.exhaustiveness),
                "--num_modes", str(self.n_poses),
                "--seed", str(self.seed),
                "--cpu", str(self.cpu),
                "--out", str(out),
                "--quiet",
                *self.extra_args(),
            ]  # fmt: skip
            run(args, cwd=work, env=self.environment())
            poses = read_poses(out)[: self.n_poses]
        return DockingResult(poses=poses, scores=[self.score_of(p) for p in poses])


class GninaEngine(SminaEngine):
    """
    Gnina: Smina's search with CNN rescoring of the poses, which also sets their rank.

    Each worker process uses one GPU, chosen from `gpus` by its process id, so parallel
    workers spread across the cards; without GPUs it runs on the CPU.

    Args:
        exhaustiveness (int): Search effort per run.
        n_poses (int): Maximum number of poses to return.
        seed (int): Random seed.
        cpu (int): Threads per docking run.
        gpus (Sequence[int] | None): GPU indices to use; None or empty runs on the CPU.
        command (Sequence[str] | None): Program to run; default looks up `gnina` on PATH.
    """

    name = "gnina"
    program = "gnina"

    def __init__(
        self,
        exhaustiveness: int = 8,
        n_poses: int = 9,
        seed: int = 42,
        cpu: int = 1,
        gpus: Sequence[int] | None = None,
        command: Sequence[str] | None = None,
    ):
        super().__init__(exhaustiveness, n_poses, seed, cpu, command)
        self.gpus = list(gpus) if gpus else []

    def device(self) -> int | None:
        """GPU for this worker process, or None for CPU."""
        return self.gpus[os.getpid() % len(self.gpus)] if self.gpus else None

    def extra_args(self) -> list[str]:
        args = ["--cnn_scoring", "rescore"]
        return args + (["--device", "0"] if self.gpus else ["--no_gpu"])

    def environment(self) -> dict | None:
        device = self.device()
        if device is None:
            return None
        # Expose only this worker's card, so the process cannot touch the others.
        return {**os.environ, "CUDA_VISIBLE_DEVICES": str(device)}
