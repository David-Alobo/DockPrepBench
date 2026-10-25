"""rDock, run as external programs (rbcavity, rbdock) with Open Babel for the receptor."""

import os
from collections.abc import Sequence
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
from rdkit import Chem
from rdkit.Geometry import Point3D

from dockprepbench.engines.base import Box, DockingEngine, DockingResult
from dockprepbench.engines.external import (
    EngineError,
    read_poses,
    resolve_command,
    run,
    write_ligand,
)

# rDock's standard docking protocol, shipped with rDock.
DOCKING_PROTOCOL = "dock.prm"


REFERENCE_FILE = "box.sdf"
MAX_REFERENCE_POINTS = 999  # V2000 SD limit; rbcavity crashes on V3000 references
MIN_SPACING = 1.5  # Å between grid points filling the box


def box_points(box: Box) -> tuple[np.ndarray, float]:
    """
    Fill a search box with a regular grid of points for rDock's ligand site mapper.

    rDock maps cavities, not boxes. Its reference-ligand mapper keeps free space within a
    radius of any reference atom, so a reference made of points filling the box, each with
    a radius that covers its grid cell, gives a cavity equal to the free space in the box.
    The grid is as fine as the 999-atom SD limit allows, and at least MIN_SPACING.

    Args:
        box (Box): Search box.

    Returns:
        tuple[np.ndarray, float]: Grid points (n, 3), one at the centre of each cell of a
            grid dividing the box, and the radius (half a cell diagonal) at which each
            point's sphere covers its cell. The cavity therefore covers the whole box and
            extends beyond a face by at most the radius minus half the cell width (about a
            third of a cell).
    """
    size = np.array(box.size)
    spacing = MIN_SPACING
    while np.prod(np.ceil(size / spacing)) > MAX_REFERENCE_POINTS:
        spacing += 0.05
    counts = np.maximum(np.ceil(size / spacing), 1).astype(int)
    cell = size / counts
    low = np.array(box.center) - size / 2
    axes = [lo + c * (np.arange(n) + 0.5) for lo, c, n in zip(low, cell, counts, strict=True)]
    points = np.array(np.meshgrid(*axes, indexing="ij")).reshape(3, -1).T
    return points, float(np.linalg.norm(cell) / 2)


def box_reference(box: Box) -> tuple[Chem.Mol, float]:
    """
    Build the box-filling reference molecule for rDock (unbonded carbon atoms).

    Args:
        box (Box): Search box.

    Returns:
        tuple[Chem.Mol, float]: Reference molecule with one conformer, and the mapper radius.
    """
    points, radius = box_points(box)
    mol = Chem.RWMol()
    for _ in points:
        mol.AddAtom(Chem.Atom(6))
    conf = Chem.Conformer(len(points))
    for i, (x, y, z) in enumerate(points):
        conf.SetAtomPosition(i, Point3D(float(x), float(y), float(z)))
    mol.AddConformer(conf, assignId=True)
    return mol.GetMol(), radius


def cavity_parameters(receptor_mol2: str, radius: float) -> str:
    """
    Write an rDock receptor parameter file whose cavity is the free space in the search box.

    The cavity is mapped by rDock's reference-ligand mapper around the box-filling points
    written to REFERENCE_FILE (see `box_reference`), so every engine searches the same box.

    Args:
        receptor_mol2 (str): Receptor file name, relative to the working directory.
        radius (float): Mapper radius from `box_reference`.

    Returns:
        str: Contents of the .prm file.
    """
    return (
        "RBT_PARAMETER_FILE_V1.00\n"
        "TITLE dockprepbench\n"
        f"RECEPTOR_FILE {receptor_mol2}\n"
        "RECEPTOR_FLEX 3.0\n"
        "SECTION MAPPER\n"
        "    SITE_MAPPER RbtLigandSiteMapper\n"
        f"    REF_MOL {REFERENCE_FILE}\n"
        f"    RADIUS {radius:.3f}\n"
        "    SMALL_SPHERE 1.0\n"
        "    MIN_VOLUME 100\n"
        "    MAX_CAVITIES 1\n"
        "    VOL_INCR 0.0\n"
        "    GRIDSTEP 0.5\n"
        "END_SECTION\n"
        "SECTION CAVITY\n"
        "    SCORING_FUNCTION RbtCavityGridSF\n"
        "    WEIGHT 1.0\n"
        "END_SECTION\n"
    )


def _graph(mol: Chem.Mol) -> Chem.Mol:
    """
    Reduce a heavy-atom molecule to elements and connectivity, for atom matching.

    Args:
        mol (Chem.Mol): Heavy-atom molecule, sanitized or not.

    Returns:
        Chem.Mol: Copy with single, non-aromatic bonds and neutral atoms; not sanitized.
    """
    graph = Chem.RWMol(mol)
    for atom in graph.GetAtoms():
        atom.SetFormalCharge(0)
        atom.SetIsAromatic(False)
        atom.SetNoImplicit(True)
        atom.SetNumExplicitHs(0)
    for bond in graph.GetBonds():
        bond.SetBondType(Chem.BondType.SINGLE)
        bond.SetIsAromatic(False)
    graph.UpdatePropertyCache(strict=False)
    return graph.GetMol()


def restore_chemistry(pose: Chem.Mol, ligand: Chem.Mol) -> Chem.Mol | None:
    """
    Put an rDock pose's heavy-atom coordinates onto the input ligand.

    rDock writes poses without formal charges, so charged groups (nitro, quaternary N) come
    back with impossible valences, and it does not keep every input hydrogen. The heavy
    atoms are matched to the input ligand on elements and connectivity alone, the input's
    heavy atoms take the pose's coordinates, and hydrogens are added back in place.

    Args:
        pose (Chem.Mol): Unsanitized pose read from rDock's output.
        ligand (Chem.Mol): The ligand as passed to rDock.

    Returns:
        Chem.Mol | None: The ligand in the pose's coordinates, with hydrogens and the pose's
            properties; None if the heavy atoms do not match the input.
    """
    heavy_pose = Chem.RemoveAllHs(pose, sanitize=False)
    template = Chem.RemoveAllHs(ligand)
    if heavy_pose.GetNumAtoms() != template.GetNumAtoms():
        return None
    match = _graph(heavy_pose).GetSubstructMatch(_graph(template))  # template atom -> pose atom
    if len(match) != template.GetNumAtoms():
        return None
    source = heavy_pose.GetConformer()
    conf = Chem.Conformer(template.GetNumAtoms())
    for i, j in enumerate(match):
        conf.SetAtomPosition(i, source.GetAtomPosition(j))
    template.RemoveAllConformers()
    template.AddConformer(conf, assignId=True)
    mol = Chem.AddHs(template, addCoords=True)
    for key, value in pose.GetPropsAsDict(includePrivate=False).items():
        mol.SetProp(key, str(value))
    return mol


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
            reference, radius = box_reference(box)
            write_ligand(reference, work / REFERENCE_FILE)
            (work / "receptor.prm").write_text(cavity_parameters("receptor.mol2", radius))
            run([*self._cmd("rbcavity"), "-W", "-r", "receptor.prm"], cwd=work, env=env)
            write_ligand(ligand, work / "ligand.sdf")
            rbdock = [
                *self._cmd("rbdock"),
                "-i", "ligand.sdf", "-o", "poses", "-r", "receptor.prm",
                "-p", DOCKING_PROTOCOL, "-n", str(self.n_runs), "-s", str(self.seed),
            ]  # fmt: skip
            run(rbdock, cwd=work, env=env)
            poses = read_poses(work / "poses.sd", sanitize=False)
        poses = [p for p in (restore_chemistry(p, ligand) for p in poses) if p is not None]
        if not poses:
            raise EngineError("rDock poses do not match the input ligand's atoms")
        poses.sort(key=lambda p: float(p.GetProp("SCORE")))
        poses = poses[: self.n_poses]
        return DockingResult(poses=poses, scores=[float(p.GetProp("SCORE")) for p in poses])
