"""AutoDock Vina through its Python bindings (Linux x86_64 only)."""

from pathlib import Path
from tempfile import TemporaryDirectory

from rdkit import Chem

from dockprepbench.engines.base import Box, DockingEngine, DockingResult
from dockprepbench.engines.pdbqt import ligand_to_pdbqt, pdbqt_to_poses, receptor_to_pdbqt


class VinaEngine(DockingEngine):
    """
    AutoDock Vina with the default Vina scoring function.

    Args:
        exhaustiveness (int): Search effort per run.
        n_poses (int): Maximum number of poses to return.
        seed (int): Random seed. Fixed by default so runs are reproducible.
        cpu (int): Number of threads; 0 uses all available.
    """

    name = "vina"

    def __init__(self, exhaustiveness: int = 8, n_poses: int = 9, seed: int = 42, cpu: int = 0):
        self.exhaustiveness = exhaustiveness
        self.n_poses = n_poses
        self.seed = seed
        self.cpu = cpu

    def dock(self, receptor: Path, ligand: Chem.Mol, box: Box) -> DockingResult:
        from vina import Vina  # deferred: the bindings only install on Linux

        vina = Vina(sf_name="vina", cpu=self.cpu, seed=self.seed, verbosity=0)
        with TemporaryDirectory() as tmp:
            receptor_pdbqt = Path(tmp) / "receptor.pdbqt"
            receptor_pdbqt.write_text(receptor_to_pdbqt(Path(receptor).read_text(), box))
            vina.set_receptor(rigid_pdbqt_filename=str(receptor_pdbqt))

        vina.set_ligand_from_string(ligand_to_pdbqt(ligand))
        vina.compute_vina_maps(center=list(box.center), box_size=list(box.size))
        vina.dock(exhaustiveness=self.exhaustiveness, n_poses=self.n_poses)

        poses = pdbqt_to_poses(vina.poses(n_poses=self.n_poses))
        scores = [float(e) for e in vina.energies(n_poses=self.n_poses)[:, 0]]
        return DockingResult(poses=poses, scores=scores)
