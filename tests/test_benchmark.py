import math
import shutil

import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from dockprepbench.benchmark import dock_complex, run_benchmark
from dockprepbench.datasets.posebusters import Complex
from dockprepbench.engines.base import DockingEngine, DockingResult


class CrystalEngine(DockingEngine):
    """Returns the crystal ligand as the top pose, so a correct pipeline scores RMSD 0."""

    name = "crystal"

    def __init__(self, crystal_path):
        self.crystal_path = crystal_path

    def dock(self, receptor, ligand, box):
        crystal = next(Chem.SDMolSupplier(str(self.crystal_path)))
        return DockingResult(poses=[crystal, ligand], scores=[-9.0, -5.0])


class EmptyEngine(DockingEngine):
    name = "empty"

    def dock(self, receptor, ligand, box):
        return DockingResult(poses=[], scores=[])


@pytest.fixture
def entry(tmp_path, peptide_pdb):
    mol = Chem.AddHs(Chem.MolFromSmiles("NC(=O)c1ccccc1"))
    crystal, start = tmp_path / "crystal.sdf", tmp_path / "start.sdf"
    for path, conf_seed in ((crystal, 1), (start, 7)):
        AllChem.EmbedMolecule(mol, randomSeed=conf_seed)
        with Chem.SDWriter(str(path)) as writer:
            writer.write(mol)
    protein = tmp_path / "protein.pdb"
    shutil.copy(peptide_pdb, protein)
    return Complex(id="TEST_LIG", protein=protein, crystal_ligand=crystal, start_ligand=start)


def test_dock_complex_scores_top_pose(entry, tmp_path):
    result = dock_complex(entry, CrystalEngine(entry.crystal_ligand), tmp_path / "work")

    assert result.top1_rmsd == pytest.approx(0.0, abs=1e-6)
    assert result.top1_score == -9.0
    assert isinstance(result.top1_pb_valid, bool)
    assert result.n_poses == 2
    assert result.error == ""
    assert (tmp_path / "work" / "TEST_LIG_receptor.pdb").is_file()


def test_dock_complex_records_failure(entry, tmp_path):
    result = dock_complex(entry, EmptyEngine(), tmp_path / "work")

    assert math.isnan(result.top1_rmsd)
    assert "no poses" in result.error


def test_run_benchmark_returns_one_row_per_complex(entry, tmp_path):
    engine = CrystalEngine(entry.crystal_ligand)

    df = run_benchmark([entry, entry], engine, tmp_path / "work")

    assert list(df["id"]) == ["TEST_LIG", "TEST_LIG"]
    assert list(df.columns) == [
        "id",
        "top1_rmsd",
        "top1_score",
        "top1_pb_valid",
        "top1_failed_checks",
        "n_poses",
        "seconds",
        "error",
    ]
