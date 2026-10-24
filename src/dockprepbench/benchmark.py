"""Run a docking engine over a set of complexes and score the top-ranked poses."""

import json
import logging
import math
import os
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict, dataclass
from pathlib import Path

import pandas as pd
from rdkit import Chem

from dockprepbench.datasets.posebusters import Complex
from dockprepbench.engines.base import DockingEngine
from dockprepbench.evaluation.rmsd import centroid_shift, pose_rmsd
from dockprepbench.evaluation.validity import is_pb_valid, validity_checks
from dockprepbench.io import load_sdf_file, write_sdf_file
from dockprepbench.preparation.box import box_from_gyration, box_from_ligand
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
    top1_centroid_shift: float = math.nan  # Å; large values mean a different site
    best_rmsd: float = math.nan  # lowest RMSD among all returned poses
    best_rank: int = 0  # 1-based rank of that pose; > 1 with best_rmsd < 2 Å = scoring miss
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
    entry: Complex,
    engine: DockingEngine,
    config: PrepConfig,
    work_dir: Path,
    pose_dir: Path | None = None,
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
        pose_dir (Path | None): If given, all returned poses are written to
            `<pose_dir>/<id>.sdf` in rank order, each with `rank`, `score` and `rmsd`
            properties, for inspection and visualisation.

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
        if config.box_size == "gyration":
            box = box_from_gyration(crystal, ligand)
        else:
            box = box_from_ligand(crystal, padding=config.box_padding)
        result = engine.dock(receptor, ligand, box)
        if not result.poses:
            raise RuntimeError(f"{engine.name} returned no poses")
        checks = validity_checks(result.poses[0], crystal, entry.protein)
        rmsds = [pose_rmsd(pose, crystal) for pose in result.poses]
        best = min(range(len(rmsds)), key=rmsds.__getitem__)
        if pose_dir is not None:
            _save_poses(result.poses, result.scores, rmsds, pose_dir / f"{entry.id}.sdf")
        return ComplexResult(
            id=entry.id,
            top1_rmsd=rmsds[0],
            top1_score=result.scores[0],
            top1_pb_valid=is_pb_valid(checks),
            top1_failed_checks=";".join(name for name, ok in checks.items() if not ok),
            top1_centroid_shift=centroid_shift(result.poses[0], crystal),
            best_rmsd=rmsds[best],
            best_rank=best + 1,
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
    pose_dir: Path | None = None,
    checkpoint: Path | None = None,
) -> pd.DataFrame:
    """
    Dock every complex and collect the results.

    With a checkpoint, each result is appended to it as soon as its complex finishes, and
    complexes already recorded there are not docked again. A run that is interrupted (or
    killed, e.g. by a memory limit) therefore resumes where it stopped when restarted.

    Args:
        complexes (list[Complex]): Complexes to dock.
        engine (DockingEngine): Engine to dock with. Must be picklable when workers > 1.
        config (PrepConfig): Preparation protocol.
        work_dir (Path): Root directory for cached prepared receptors.
        workers (int): Number of complexes docked in parallel.
        pose_dir (Path | None): If given, every complex's poses are written there as
            `<id>.sdf` (see `dock_complex`).
        checkpoint (Path | None): JSON Lines file of finished results, one per line.

    Returns:
        pd.DataFrame: One row per complex, in input order, with the ComplexResult fields.
    """
    wanted = {c.id for c in complexes}
    done = {
        k: v for k, v in (_load_checkpoint(checkpoint) if checkpoint else {}).items() if k in wanted
    }
    todo = [c for c in complexes if c.id not in done]
    total = len(complexes)
    if done:
        logger.info("Resuming: %d of %d complexes already done", len(done), total)

    def record(result: ComplexResult) -> None:
        done[result.id] = asdict(result)
        if checkpoint:
            _append_checkpoint(checkpoint, result)
        if len(done) % 10 == 0 or len(done) == total:
            logger.info("Progress: %d of %d complexes done", len(done), total)

    args = (engine, config, work_dir, pose_dir)
    if workers == 1:
        for entry in todo:
            record(dock_complex(entry, *args))
    else:
        with ProcessPoolExecutor(workers) as pool:
            futures = [pool.submit(dock_complex, entry, *args) for entry in todo]
            for future in as_completed(futures):
                record(future.result())
    return pd.DataFrame([done[c.id] for c in complexes])


def _load_checkpoint(path: Path) -> dict[str, dict]:
    """
    Read finished results from a checkpoint, skipping a final line cut off by a crash.

    Args:
        path (Path): JSON Lines checkpoint; may not exist yet.

    Returns:
        dict[str, dict]: ComplexResult fields keyed by complex id.
    """
    if not path.is_file():
        return {}
    done = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        try:
            row = json.loads(line)
        except json.JSONDecodeError:
            logger.warning("Ignoring an incomplete line in %s", path)
            continue
        done[row["id"]] = row
    return done


def _append_checkpoint(path: Path, result: ComplexResult) -> None:
    """
    Append one result to the checkpoint and flush it to disk immediately.

    Args:
        path (Path): JSON Lines checkpoint.
        result (ComplexResult): Finished result.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as f:
        f.write(json.dumps(asdict(result)) + "\n")
        f.flush()
        os.fsync(f.fileno())


def _save_poses(poses: list[Chem.Mol], scores: list[float], rmsds: list[float], path: Path) -> None:
    """
    Write ranked poses to an SDF file, tagging each with its rank, score and RMSD.

    Args:
        poses (list[Chem.Mol]): Poses in rank order.
        scores (list[float]): Engine score per pose.
        rmsds (list[float]): RMSD to the crystal ligand per pose, in Å.
        path (Path): Destination SDF file.
    """
    tagged = []
    for rank, (pose, score, rmsd) in enumerate(zip(poses, scores, rmsds, strict=True), 1):
        mol = Chem.Mol(pose)
        mol.SetProp("_Name", f"{path.stem} rank {rank}")
        mol.SetIntProp("rank", rank)
        mol.SetDoubleProp("score", score)
        mol.SetDoubleProp("rmsd", rmsd)
        tagged.append(mol)
    write_sdf_file(tagged, path)
