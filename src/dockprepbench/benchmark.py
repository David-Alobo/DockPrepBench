"""Run a docking engine over a set of complexes and score the top-ranked poses."""

import logging
import math
import time
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict, dataclass
from pathlib import Path

import pandas as pd

from dockprepbench.datasets.posebusters import Complex
from dockprepbench.engines.base import DockingEngine
from dockprepbench.evaluation.rmsd import pose_rmsd
from dockprepbench.evaluation.validity import is_pb_valid, validity_checks
from dockprepbench.io import load_sdf_file
from dockprepbench.preparation.box import box_from_ligand
from dockprepbench.preparation.config import PrepConfig
from dockprepbench.preparation.ligand import prepare_ligand
from dockprepbench.preparation.receptor import protonate_receptor

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class ComplexResult:
    """Outcome of docking one complex. Failed runs have NaN metrics and an error message."""

    id: str
    top1_rmsd: float = math.nan
    top1_score: float = math.nan
    top1_pb_valid: bool = False
    top1_failed_checks: str = ""  # semicolon-separated PoseBusters check names
    n_poses: int = 0
    seconds: float = math.nan
    error: str = ""


def prepared_receptor(entry: Complex, config: PrepConfig, work_dir: Path) -> Path:
    """
    Return the protonated receptor for a complex, preparing and caching it on first use.

    Receptors are cached per `config.receptor_key()`, so configs that differ only in
    settings that do not touch the receptor share them.

    Args:
        entry (Complex): Complex whose protein to prepare.
        config (PrepConfig): Preparation protocol.
        work_dir (Path): Root directory for cached prepared receptors.

    Returns:
        Path: Prepared receptor PDB file.
    """
    cache_dir = work_dir / config.receptor_key()
    path = cache_dir / f"{entry.id}_receptor.pdb"
    if not path.is_file():
        cache_dir.mkdir(parents=True, exist_ok=True)
        tmp = path.with_suffix(".tmp")
        tmp.write_text(protonate_receptor(entry.protein, ph=config.ph))
        tmp.replace(path)  # atomic, so an interrupted run never leaves a partial file
    return path


def dock_complex(
    entry: Complex, engine: DockingEngine, config: PrepConfig, work_dir: Path
) -> ComplexResult:
    """
    Prepare and redock a complex's ligand, starting from its supplied start conformer,
    and score the top pose.

    The top pose gets an RMSD against the crystal ligand and the PoseBusters validity
    checks, run against the original protein file (with cofactors), as in PoseBusters.

    Any error is caught and recorded so one bad complex does not stop a benchmark run.

    Args:
        entry (Complex): Complex to dock.
        engine (DockingEngine): Engine to dock with.
        config (PrepConfig): Preparation protocol.
        work_dir (Path): Root directory for cached prepared receptors.

    Returns:
        ComplexResult: Top-1 RMSD against the crystal ligand, score and timing.
    """
    start = time.perf_counter()
    try:
        crystal = load_sdf_file(entry.crystal_ligand)[0]
        ligand = prepare_ligand(
            load_sdf_file(entry.start_ligand)[0],
            protonation=config.ligand_protonation,
            conformer=config.ligand_conformer,
            ph=config.ph,
        )
        receptor = prepared_receptor(entry, config, work_dir)
        box = box_from_ligand(crystal, padding=config.box_padding)
        result = engine.dock(receptor, ligand, box)
        if not result.poses:
            raise RuntimeError(f"{engine.name} returned no poses")
        checks = validity_checks(result.poses[0], crystal, entry.protein)
        return ComplexResult(
            id=entry.id,
            top1_rmsd=pose_rmsd(result.poses[0], crystal),
            top1_score=result.scores[0],
            top1_pb_valid=is_pb_valid(checks),
            top1_failed_checks=";".join(name for name, ok in checks.items() if not ok),
            n_poses=len(result.poses),
            seconds=time.perf_counter() - start,
        )
    except Exception as err:
        logger.warning("%s failed: %s: %s", entry.id, type(err).__name__, err)
        return ComplexResult(
            id=entry.id, seconds=time.perf_counter() - start, error=f"{type(err).__name__}: {err}"
        )


def run_benchmark(
    complexes: list[Complex],
    engine: DockingEngine,
    config: PrepConfig,
    work_dir: Path,
    workers: int = 1,
) -> pd.DataFrame:
    """
    Dock every complex and collect the results.

    Args:
        complexes (list[Complex]): Complexes to dock.
        engine (DockingEngine): Engine to dock with. Must be picklable when workers > 1.
        config (PrepConfig): Preparation protocol.
        work_dir (Path): Root directory for cached prepared receptors.
        workers (int): Number of complexes docked in parallel.

    Returns:
        pd.DataFrame: One row per complex, in input order, with the ComplexResult fields.
    """
    n = len(complexes)
    args = (complexes, [engine] * n, [config] * n, [work_dir] * n)
    if workers == 1:
        results = list(map(dock_complex, *args))
    else:
        with ProcessPoolExecutor(workers) as pool:
            results = list(pool.map(dock_complex, *args))
    return pd.DataFrame([asdict(r) for r in results])
