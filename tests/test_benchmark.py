import json
import math
import shutil
from dataclasses import asdict

import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from dockprepbench.benchmark import ComplexResult, dock_complex, run_benchmark
from dockprepbench.datasets.posebusters import Complex
from dockprepbench.engines.base import DockingEngine, DockingResult
from dockprepbench.preparation.config import PrepConfig


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
    result = dock_complex(
        entry, CrystalEngine(entry.crystal_ligand), PrepConfig(), tmp_path / "work"
    )

    assert result.top1_rmsd == pytest.approx(0.0, abs=1e-6)
    assert result.top1_score == -9.0
    assert isinstance(result.top1_pb_valid, bool)
    assert result.top1_centroid_shift == pytest.approx(0.0, abs=1e-6)
    assert (result.best_rmsd, result.best_rank) == (pytest.approx(0.0, abs=1e-6), 1)
    assert result.n_poses == 2
    assert result.error == ""
    assert (tmp_path / "work" / "ph7.4" / "TEST_LIG_receptor.pdb").is_file()


def test_dock_complex_records_failure(entry, tmp_path):
    result = dock_complex(entry, EmptyEngine(), PrepConfig(), tmp_path / "work")

    assert math.isnan(result.top1_rmsd)
    assert "no poses" in result.error


def test_run_benchmark_returns_one_row_per_complex(entry, tmp_path):
    engine = CrystalEngine(entry.crystal_ligand)

    df = run_benchmark([entry, entry], engine, PrepConfig(), tmp_path / "work")

    assert list(df["id"]) == ["TEST_LIG", "TEST_LIG"]
    assert list(df.columns) == [
        "id",
        "top1_rmsd",
        "top1_score",
        "top1_pb_valid",
        "top1_failed_checks",
        "top1_centroid_shift",
        "best_rmsd",
        "best_rank",
        "n_poses",
        "seconds",
        "error",
    ]


class ReversedEngine(CrystalEngine):
    """Ranks the start conformer first and the crystal pose second: a scoring miss."""

    def dock(self, receptor, ligand, box):
        result = super().dock(receptor, ligand, box)
        return DockingResult(poses=result.poses[::-1], scores=result.scores[::-1])


def test_dock_complex_records_best_pose_below_top(entry, tmp_path):
    result = dock_complex(
        entry, ReversedEngine(entry.crystal_ligand), PrepConfig(), tmp_path / "work"
    )

    assert result.top1_rmsd > 0.1
    assert result.best_rmsd == pytest.approx(0.0, abs=1e-6)
    assert result.best_rank == 2


def test_dock_complex_saves_ranked_poses(entry, tmp_path):
    pose_dir = tmp_path / "poses"

    dock_complex(
        entry, ReversedEngine(entry.crystal_ligand), PrepConfig(), tmp_path / "work", pose_dir
    )

    poses = list(Chem.SDMolSupplier(str(pose_dir / "TEST_LIG.sdf")))
    assert [p.GetIntProp("rank") for p in poses] == [1, 2]
    assert [p.GetDoubleProp("score") for p in poses] == [-5.0, -9.0]
    assert poses[1].GetDoubleProp("rmsd") == pytest.approx(0.0, abs=1e-6)


def test_failed_complex_writes_no_poses(entry, tmp_path):
    pose_dir = tmp_path / "poses"

    dock_complex(entry, EmptyEngine(), PrepConfig(), tmp_path / "work", pose_dir)

    assert not pose_dir.exists()


def test_run_benchmark_checkpoints_every_result(entry, tmp_path):
    checkpoint = tmp_path / "out" / "results.partial.jsonl"
    other = Complex(
        id="OTHER_LIG",
        protein=entry.protein,
        crystal_ligand=entry.crystal_ligand,
        start_ligand=entry.start_ligand,
    )

    df = run_benchmark(
        [entry, other],
        CrystalEngine(entry.crystal_ligand),
        PrepConfig(),
        tmp_path / "work",
        checkpoint=checkpoint,
    )

    lines = checkpoint.read_text(encoding="utf-8").splitlines()
    assert len(lines) == 2
    assert list(df["id"]) == ["TEST_LIG", "OTHER_LIG"]


def test_run_benchmark_resumes_without_redocking(entry, tmp_path):
    checkpoint = tmp_path / "results.partial.jsonl"
    finished = ComplexResult(id="TEST_LIG", top1_rmsd=0.5, top1_score=-9.0, n_poses=9)
    checkpoint.write_text(json.dumps(asdict(finished)) + "\n", encoding="utf-8")

    # EmptyEngine would record an error if TEST_LIG were docked again.
    df = run_benchmark(
        [entry], EmptyEngine(), PrepConfig(), tmp_path / "work", checkpoint=checkpoint
    )

    assert df.loc[0, "top1_rmsd"] == 0.5
    assert df.loc[0, "error"] == ""


def test_run_benchmark_ignores_a_line_cut_off_by_a_crash(entry, tmp_path):
    checkpoint = tmp_path / "results.partial.jsonl"
    finished = ComplexResult(id="TEST_LIG", top1_rmsd=0.5, n_poses=9)
    checkpoint.write_text(
        json.dumps(asdict(finished)) + "\n" + '{"id": "OTHER_LIG", "top1_rm', encoding="utf-8"
    )

    df = run_benchmark(
        [entry], EmptyEngine(), PrepConfig(), tmp_path / "work", checkpoint=checkpoint
    )

    assert df.loc[0, "top1_rmsd"] == 0.5


def test_run_benchmark_ignores_checkpointed_complexes_not_requested(entry, tmp_path):
    checkpoint = tmp_path / "results.partial.jsonl"
    stray = ComplexResult(id="NOT_IN_RUN", top1_rmsd=0.5)
    checkpoint.write_text(json.dumps(asdict(stray)) + "\n", encoding="utf-8")

    df = run_benchmark(
        [entry],
        CrystalEngine(entry.crystal_ligand),
        PrepConfig(),
        tmp_path / "work",
        checkpoint=checkpoint,
    )

    assert list(df["id"]) == ["TEST_LIG"]
