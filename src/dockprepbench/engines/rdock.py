"""rDock, run as external programs (rbcavity, rbdock) with Open Babel for the receptor."""

import os
from collections.abc import Sequence
from pathlib import Path
from tempfile import TemporaryDirectory

from rdkit import Chem

from dockprepbench.engines.base import Box, DockingEngine, DockingResult
from dockprepbench.engines.external import read_poses, resolve_command, run, write_ligand

# rDock's standard docking protocol, shipped with rDock.
DOCKING_PROTOCOL = "dock.prm"


def cavity_parameters(receptor_mol2: str, box: Box) -> str:
    """
    Write an rDock receptor parameter file whose cavity is a sphere enclosing the box.

    rDock maps cavities rather than boxes; a sphere through the box's corners keeps the
    whole search box available, so the box-padding arms mean the same thing as for the
    other engines.

    Args:
        receptor_mol2 (str): Receptor file name, relative to the working directory.
        box (Box): Search box.

    Returns:
        str: Contents of the .prm file.
    """
    cx, cy, cz = box.center
    radius = 0.5 * sum(s * s for s in box.size) ** 0.5
    return (
        "RBT_PARAMETER_FILE_V1.00\n"
        "TITLE dockprepbench\n"
        f"RECEPTOR_FILE {receptor_mol2}\n"
        "RECEPTOR_FLEX 3.0\n"
        "SECTION MAPPER\n"
        "    SITE_MAPPER RbtSphereSiteMapper\n"
        f"    CENTER ({cx:.3f},{cy:.3f},{cz:.3f})\n"
        f"    RADIUS {radius:.3f}\n"
        "    SMALL_SPHERE 1.0\n"
        "    LARGE_SPHERE 3.5\n"
        "    MAX_CAVITIES 1\n"
        "END_SECTION\n"
        "SECTION CAVITY\n"
        "    SCORING_FUNCTION RbtCavityGridSF\n"
        "    WEIGHT 1.0\n"
        "END_SECTION\n"
    )


class RDockEngine(DockingEngine):
    """
    rDock with its standard docking protocol; poses are ranked by total SCORE.

    Args:
        n_runs (int): Independent docking runs per ligand (rDock's recommendation is 50).
        n_poses (int): Number of best-scoring poses to return.
        seed (int): Random seed.
        commands (dict[str, Sequence[str]] | None): Overrides for `rbcavity`, `rbdock` and
            `obabel`; defaults look them up on PATH.
    """

    name = "rdock"

    def __init__(
        self,
        n_runs: int = 50,
        n_poses: int = 9,
        seed: int = 42,
        commands: dict[str, Sequence[str]] | None = None,
    ):
        self.n_runs = n_runs
        self.n_poses = n_poses
        self.seed = seed
        self.commands = {k: list(v) for k, v in (commands or {}).items()}

    def _cmd(self, program: str) -> list[str]:
        """
        Command prefix for one of rDock's programs (or Open Babel).

        Args:
            program (str): `rbcavity`, `rbdock` or `obabel`.

        Returns:
            list[str]: Override from `commands`, or the program found on PATH.
        """
        return resolve_command(self.commands.get(program), program)

    def _environment(self) -> dict:
        """
        Environment for rDock: RBT_ROOT must point at rDock's data directory.

        When unset it is inferred from where `rbdock` lives: a release tarball
        (`<root>/bin/rbdock` beside `<root>/data`, whose `<root>/lib` is also added to the
        library path) or a conda environment (`<env>/share/rdock-*`).

        Returns:
            dict: Environment for rDock's processes.
        """
        env = dict(os.environ)
        if "RBT_ROOT" in env:
            return env
        prefix = Path(self._cmd("rbdock")[-1]).resolve().parent.parent
        if (prefix / "data" / "scripts" / DOCKING_PROTOCOL).is_file():
            env["RBT_ROOT"] = str(prefix)
            lib = str(prefix / "lib")
            env["LD_LIBRARY_PATH"] = os.pathsep.join(
                filter(None, [lib, env.get("LD_LIBRARY_PATH")])
            )
        elif roots := sorted((prefix / "share").glob("rdock*")):
            env["RBT_ROOT"] = str(roots[-1])
        return env

    def dock(self, receptor: Path, ligand: Chem.Mol, box: Box) -> DockingResult:
        env = self._environment()
        with TemporaryDirectory() as tmp:
            work = Path(tmp)
            run(
                [*self._cmd("obabel"), str(Path(receptor).resolve()), "-O", "receptor.mol2"],
                cwd=work,
            )
            (work / "receptor.prm").write_text(cavity_parameters("receptor.mol2", box))
            run([*self._cmd("rbcavity"), "-W", "-r", "receptor.prm"], cwd=work, env=env)
            write_ligand(ligand, work / "ligand.sdf")
            rbdock = [
                *self._cmd("rbdock"),
                "-i", "ligand.sdf", "-o", "poses", "-r", "receptor.prm",
                "-p", DOCKING_PROTOCOL, "-n", str(self.n_runs), "-s", str(self.seed),
            ]  # fmt: skip
            run(rbdock, cwd=work, env=env)
            poses = read_poses(work / "poses.sd")
        poses.sort(key=lambda p: float(p.GetProp("SCORE")))
        poses = poses[: self.n_poses]
        return DockingResult(poses=poses, scores=[float(p.GetProp("SCORE")) for p in poses])
