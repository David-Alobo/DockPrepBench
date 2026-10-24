import json
import shutil
import sys
from pathlib import Path

import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from dockprepbench.engines.base import Box
from dockprepbench.engines.external import EngineError, resolve_command
from dockprepbench.engines.rdock import RDockEngine, cavity_parameters
from dockprepbench.engines.smina import GninaEngine, SminaEngine

FAKE = [sys.executable, str(Path(__file__).with_name("fake_engine.py"))]
BOX = Box(center=(1.0, 2.0, 3.0), size=(10.0, 12.0, 14.0))


@pytest.fixture
def ligand():
    mol = Chem.AddHs(Chem.MolFromSmiles("NC(=O)c1ccccc1"))
    AllChem.EmbedMolecule(mol, randomSeed=1)
    return mol


@pytest.fixture
def calls(tmp_path, monkeypatch):
    """Path where the fake engine copies the command lines it received."""
    log = tmp_path / "calls.json"
    monkeypatch.setenv("FAKE_ENGINE_LOG", str(log))
    monkeypatch.delenv("FAKE_ENGINE_FAIL", raising=False)
    return lambda: json.loads(log.read_text())


def arg(argv, flag):
    return argv[argv.index(flag) + 1]


def test_smina_passes_box_and_settings(peptide_pdb, ligand, calls):
    result = SminaEngine(exhaustiveness=4, n_poses=3, seed=7, command=FAKE).dock(
        peptide_pdb, ligand, BOX
    )

    (call,) = calls()
    argv = call["args"]
    assert (arg(argv, "--center_x"), arg(argv, "--size_z")) == ("1.000", "14.000")
    assert (arg(argv, "--exhaustiveness"), arg(argv, "--num_modes"), arg(argv, "--seed")) == (
        "4",
        "3",
        "7",
    )
    assert "--cnn_scoring" not in argv
    assert result.scores == [-9.0, -8.0, -7.0]
    assert Chem.MolToSmiles(Chem.RemoveHs(result.poses[0])) == "NC(=O)c1ccccc1"


def test_gnina_rescores_with_cnn_on_its_assigned_gpu(peptide_pdb, ligand, calls):
    engine = GninaEngine(n_poses=2, gpus=[3, 5], command=FAKE)

    engine.dock(peptide_pdb, ligand, BOX)

    (call,) = calls()
    assert arg(call["args"], "--cnn_scoring") == "rescore"
    assert arg(call["args"], "--device") == "0"  # the only card the process can see
    assert call["cuda"] in {"3", "5"}
    assert call["cuda"] == str(engine.device())


def test_gnina_without_gpus_runs_on_cpu(peptide_pdb, ligand, calls):
    GninaEngine(n_poses=1, command=FAKE).dock(peptide_pdb, ligand, BOX)

    (call,) = calls()
    assert "--no_gpu" in call["args"]
    assert "--device" not in call["args"]


def test_engine_failure_is_reported(peptide_pdb, ligand, calls, monkeypatch):
    monkeypatch.setenv("FAKE_ENGINE_FAIL", "1")

    with pytest.raises(EngineError, match="exited 3: simulated failure"):
        SminaEngine(command=FAKE).dock(peptide_pdb, ligand, BOX)


def test_missing_program_is_reported():
    with pytest.raises(EngineError, match="not found on PATH"):
        resolve_command(None, "definitely-not-a-docking-program")


def test_rdock_runs_cavity_then_docking_and_sorts_by_score(peptide_pdb, ligand, calls):
    engine = RDockEngine(
        n_runs=5, n_poses=2, seed=11, commands={p: FAKE for p in ("obabel", "rbcavity", "rbdock")}
    )

    result = engine.dock(peptide_pdb, ligand, BOX)

    obabel, rbcavity, rbdock = (c["args"] for c in calls())
    assert obabel[-2:] == ["-O", "receptor.mol2"]
    assert "receptor.prm" in rbcavity
    assert (arg(rbdock, "-n"), arg(rbdock, "-s"), arg(rbdock, "-p")) == ("5", "11", "dock.prm")
    assert result.scores == [-12.0, -8.0]  # best first, cut to n_poses


def test_rdock_cavity_sphere_encloses_the_box():
    prm = cavity_parameters("receptor.mol2", BOX)

    assert "CENTER (1.000,2.000,3.000)" in prm
    radius = float(next(line.split()[1] for line in prm.splitlines() if "RADIUS" in line))
    assert radius == pytest.approx(0.5 * (10**2 + 12**2 + 14**2) ** 0.5, abs=1e-3)


real_engines = {
    "smina": (SminaEngine, ["smina"]),
    "gnina": (GninaEngine, ["gnina"]),
}


@pytest.mark.parametrize("name", list(real_engines))
def test_real_engine_docks_into_box(name, peptide_pdb, ligand):
    cls, programs = real_engines[name]
    if not all(shutil.which(p) for p in programs):
        pytest.skip(f"{name} not installed")
    from dockprepbench.io import load_pdb_file

    atoms = list(load_pdb_file(peptide_pdb).get_atoms())
    centre = tuple(sum(a.coord[i] for a in atoms) / len(atoms) for i in range(3))
    engine = cls(n_runs=3) if name == "rdock" else cls(exhaustiveness=1)

    result = engine.dock(peptide_pdb, ligand, Box(center=centre, size=(20.0, 20.0, 20.0)))

    assert 1 <= len(result.poses) <= 9
    assert result.scores == sorted(result.scores) or name == "gnina"  # gnina ranks by CNN
    for pose in result.poses:
        assert Chem.MolToSmiles(Chem.RemoveHs(pose)) == "NC(=O)c1ccccc1"


POSEBUSTERS = (
    Path(__file__).parents[2] / "data" / "raw" / "posebusters" / "posebusters_benchmark_set"
)


def test_real_rdock_redocks_a_posebusters_complex(tmp_path):
    # rDock maps buried cavities, so it needs a real pocket rather than the toy peptide.
    if not all(shutil.which(p) for p in ("rbdock", "rbcavity", "obabel")):
        pytest.skip("rdock not installed")
    entry = POSEBUSTERS / "5SAK_ZRY"
    if not entry.is_dir():
        pytest.skip("PoseBusters data not downloaded")
    from dockprepbench.evaluation.rmsd import pose_rmsd
    from dockprepbench.io import load_sdf_file
    from dockprepbench.preparation.box import box_from_ligand
    from dockprepbench.preparation.receptor import protonate_receptor

    receptor = tmp_path / "receptor.pdb"
    receptor.write_text(protonate_receptor(entry / "5SAK_ZRY_protein.pdb"))
    crystal = load_sdf_file(entry / "5SAK_ZRY_ligand.sdf")[0]
    start = load_sdf_file(entry / "5SAK_ZRY_ligand_start_conf.sdf")[0]

    result = RDockEngine(n_runs=10).dock(receptor, start, box_from_ligand(crystal))

    assert 1 <= len(result.poses) <= 9
    assert result.scores == sorted(result.scores)
    assert all(pose_rmsd(p, crystal) < 30 for p in result.poses)  # poses stay in the site
