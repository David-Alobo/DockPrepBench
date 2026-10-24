"""Command-line interface: `dockprepbench`."""

from __future__ import annotations

import argparse
import json
import logging
from datetime import UTC, datetime
from pathlib import Path
from typing import TYPE_CHECKING

from dockprepbench import __version__
from dockprepbench.log import setup_logging

if TYPE_CHECKING:  # heavy imports stay lazy so `--help` and `summarize` start fast
    from dockprepbench.datasets.posebusters import Complex
    from dockprepbench.engines.base import DockingEngine
    from dockprepbench.preparation.config import PrepConfig

logger = logging.getLogger(__name__)

RESULTS_FILE = "results.csv"
CHECKPOINT_FILE = "results.partial.jsonl"  # per-complex results of an unfinished run
DEFAULT_SEED = 42  # matches VinaEngine's default, so unseeded runs stay comparable


def build_parser() -> argparse.ArgumentParser:
    """
    Build the `dockprepbench` argument parser.

    Returns:
        argparse.ArgumentParser: Parser with the global options and subcommands.
    """
    parser = argparse.ArgumentParser(
        prog="dockprepbench",
        description="Benchmark how preparation choices affect classical docking.",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable debug logging.")
    commands = parser.add_subparsers(dest="command")

    run = commands.add_parser("run", help="Redock the PoseBusters set with one engine.")
    run.add_argument("--engine", choices=["vina", "smina", "gnina", "rdock"], default="vina")
    run.add_argument(
        "--gpus",
        default="",
        help="Comma-separated GPU indices for gnina, spread across workers (e.g. 0,1,2). "
        "Empty runs gnina on the CPU.",
    )
    run.add_argument("--runs", type=int, default=50, help="rDock docking runs per ligand.")
    run.add_argument(
        "--config",
        type=Path,
        nargs="+",
        help="Preparation config(s) (TOML), run one after another. Defaults to the baseline.",
    )
    run.add_argument(
        "--force", action="store_true", help="Rerun configs whose results already exist."
    )
    run.add_argument("--data-dir", type=Path, default=Path("data/raw"))
    run.add_argument("--work-dir", type=Path, default=Path("data/prepared"))
    run.add_argument(
        "--out",
        type=Path,
        help="Run directory, for a single config. Defaults to results/<engine>/<config name>.",
    )
    run.add_argument("--limit", type=int, help="Dock only the first N complexes.")
    run.add_argument("--workers", type=int, default=1, help="Complexes docked in parallel.")
    run.add_argument("--exhaustiveness", type=int, default=8)
    run.add_argument(
        "--seed",
        type=int,
        help=f"Docking seed (default {DEFAULT_SEED}). An explicit seed adds _seed<N> to the "
        "default run directory, so repeat runs for measuring noise never overwrite each other.",
    )

    summarize = commands.add_parser("summarize", help="Success rates of one run.")
    summarize.add_argument("run", type=Path, help="Run directory or results CSV.")

    compare = commands.add_parser(
        "compare", help="Paired comparison of one or more runs against a reference run."
    )
    compare.add_argument("run_a", type=Path, help="Reference run directory or results CSV.")
    compare.add_argument(
        "run_b", type=Path, nargs="+", help="Run(s) to compare; p-values are Holm-adjusted."
    )
    compare.add_argument(
        "--pb-valid", action="store_true", help="Require PB-valid poses for success."
    )

    report = commands.add_parser(
        "report", help="Regenerate every results table (all engines and arms) as CSV files."
    )
    report.add_argument("--results", type=Path, default=Path("results"), help="Results root.")
    report.add_argument("--out", type=Path, default=Path("paper/tables"), help="Output directory.")
    report.add_argument(
        "--dataset",
        type=Path,
        default=Path("data/raw/posebusters"),
        help="Extracted PoseBusters set, for the sensitivity tables; skipped if absent.",
    )

    figures = commands.add_parser("figures", help="Figures and summary table for runs.")
    figures.add_argument(
        "runs", type=Path, nargs="+", help="Run directories or CSVs; the first is the reference."
    )
    figures.add_argument("--out", type=Path, default=Path("figures"), help="Output directory.")
    return parser


def run_command(args: argparse.Namespace) -> int:
    """
    Execute the `run` subcommand for one or more preparation configs.

    Every config is loaded before any docking starts, so a bad file fails fast. Configs
    whose results already exist are skipped unless `--force` is given, so an
    interrupted batch can simply be restarted.

    Args:
        args (argparse.Namespace): Parsed `run` arguments.

    Returns:
        int: Process exit code.
    """
    from dockprepbench.datasets.posebusters import download, load_complexes
    from dockprepbench.preparation.config import PrepConfig, load_config

    configs = [load_config(path) for path in args.config] if args.config else [PrepConfig()]
    if args.out and len(configs) > 1:
        raise SystemExit("--out takes a single config; omit it to use results/<engine>/<name>")
    names = [config.name for config in configs]
    if len(set(names)) != len(names):
        raise SystemExit(f"Config names must be unique, got {names}")

    complexes = load_complexes(download(args.data_dir))[: args.limit]
    # One Vina thread per complex when complexes already run in parallel.
    seed = DEFAULT_SEED if args.seed is None else args.seed
    suffix = "" if args.seed is None else f"_seed{args.seed}"
    engine = _make_engine(args, seed)
    for i, config in enumerate(configs, 1):
        out = args.out or Path("results") / args.engine / f"{config.name}{suffix}"
        if (out / RESULTS_FILE).is_file() and not args.force:
            logger.info("[%d/%d] %s: results exist, skipping", i, len(configs), config.name)
            continue
        logger.info(
            "[%d/%d] %s: docking %d complexes", i, len(configs), config.name, len(complexes)
        )
        _run_one(complexes, engine, config, out, args)
        print(_format_summary(out / RESULTS_FILE))
    return 0


def _make_engine(args: argparse.Namespace, seed: int) -> DockingEngine:
    """
    Build the docking engine selected on the command line.

    With several workers each docking run gets one thread, since complexes already run
    in parallel.

    Args:
        args (argparse.Namespace): Parsed `run` arguments.
        seed (int): Docking seed.

    Returns:
        DockingEngine: The configured engine.
    """
    threads = 1 if args.workers > 1 else 0
    if args.engine == "vina":
        from dockprepbench.engines.vina import VinaEngine

        return VinaEngine(exhaustiveness=args.exhaustiveness, seed=seed, cpu=threads)
    if args.engine == "smina":
        from dockprepbench.engines.smina import SminaEngine

        return SminaEngine(exhaustiveness=args.exhaustiveness, seed=seed, cpu=max(threads, 1))
    if args.engine == "gnina":
        from dockprepbench.engines.smina import GninaEngine

        gpus = [int(g) for g in args.gpus.split(",") if g.strip()]
        return GninaEngine(
            exhaustiveness=args.exhaustiveness, seed=seed, cpu=max(threads, 1), gpus=gpus
        )
    from dockprepbench.engines.rdock import RDockEngine

    return RDockEngine(n_runs=args.runs, seed=seed)


def _run_one(
    complexes: list[Complex],
    engine: DockingEngine,
    config: PrepConfig,
    out: Path,
    args: argparse.Namespace,
) -> None:
    """
    Dock all complexes under one config; write `poses/`, `run.json` and `results.csv`.

    Args:
        complexes (list[Complex]): Complexes to dock.
        engine (DockingEngine): Engine to dock with.
        config (PrepConfig): Preparation protocol.
        out (Path): Run directory.
        args (argparse.Namespace): Parsed `run` arguments (work dir, workers).
    """
    from dockprepbench.benchmark import run_benchmark

    started = datetime.now(UTC)
    out.mkdir(parents=True, exist_ok=True)
    checkpoint = out / CHECKPOINT_FILE
    results = run_benchmark(
        complexes,
        engine,
        config,
        args.work_dir,
        workers=args.workers,
        pose_dir=out / "poses",
        checkpoint=checkpoint,
    )
    metadata = {
        "dockprepbench_version": __version__,
        "engine": engine.name,
        "engine_settings": vars(engine),
        "config": config.to_dict(),
        "n_complexes": len(complexes),
        "started": started.isoformat(timespec="seconds"),
        "finished": datetime.now(UTC).isoformat(timespec="seconds"),
    }
    (out / "run.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    # Written last: its presence marks the run as complete for the skip check.
    results.to_csv(out / RESULTS_FILE, index=False)
    checkpoint.unlink(missing_ok=True)
    logger.info("Results written to %s", out)


def summarize_command(args: argparse.Namespace) -> int:
    """
    Execute the `summarize` subcommand.

    Args:
        args (argparse.Namespace): Parsed `summarize` arguments.

    Returns:
        int: Process exit code.
    """
    print(_format_summary(_results_path(args.run)))
    return 0


def compare_command(args: argparse.Namespace) -> int:
    """
    Execute the `compare` subcommand: each run B against the reference run A.

    With several B runs, p-values are also Holm-adjusted across the comparisons, since
    testing many arms against one baseline would otherwise yield false positives.

    Args:
        args (argparse.Namespace): Parsed `compare` arguments.

    Returns:
        int: Process exit code.
    """
    from dockprepbench.analysis.stats import (
        compare,
        holm_adjust,
        load_run,
        success_rate,
        successes,
    )

    path_a = _results_path(args.run_a)
    hits_a = successes(load_run(path_a), pb_valid=args.pb_valid)
    comparisons = []
    for run_b in args.run_b:
        path_b = _results_path(run_b)
        hits_b = successes(load_run(path_b), pb_valid=args.pb_valid)
        comparisons.append((path_b, hits_b, compare(hits_a, hits_b)))
    adjusted = holm_adjust([result.p_value for _, _, result in comparisons])

    metric = "PB-valid success" if args.pb_valid else "Top-1 success"
    blocks = []
    for (path_b, hits_b, result), p_holm in zip(comparisons, adjusted, strict=True):
        shared = hits_a.index.intersection(hits_b.index)
        lines = [
            f"{metric} on {result.n} shared complexes",
            f"  A {path_a}: {success_rate(hits_a[shared])}",
            f"  B {path_b}: {success_rate(hits_b[shared])}",
            f"  both {result.both}, neither {result.neither}, "
            f"only A {len(result.only_a)}, only B {len(result.only_b)}",
            f"  McNemar exact p = {result.p_value:.3g}",
        ]
        if len(comparisons) > 1:
            lines.append(f"  Holm-adjusted p = {p_holm:.3g} ({len(comparisons)} comparisons)")
        if result.only_a:
            lines.append(f"  only A: {' '.join(result.only_a)}")
        if result.only_b:
            lines.append(f"  only B: {' '.join(result.only_b)}")
        blocks.append("\n".join(lines))
    print("\n\n".join(blocks))
    return 0


def report_command(args: argparse.Namespace) -> int:
    """
    Execute the `report` subcommand: write every results table as a CSV file.

    Args:
        args (argparse.Namespace): Parsed `report` arguments.

    Returns:
        int: Process exit code.
    """
    from dockprepbench.analysis.report import write_report
    from dockprepbench.datasets.posebusters import SET_NAME

    dataset = args.dataset if (args.dataset / SET_NAME).is_dir() else None
    if dataset is None:
        print(f"No dataset at {args.dataset}; skipping the sensitivity tables.")
    for path in write_report(args.results, args.out, dataset=dataset):
        print(path)
    return 0


def figures_command(args: argparse.Namespace) -> int:
    """
    Execute the `figures` subcommand: success-rate forest plot, RMSD curves and table.

    Args:
        args (argparse.Namespace): Parsed `figures` arguments.

    Returns:
        int: Process exit code.
    """
    from dockprepbench.analysis.figures import plot_rmsd_cdf, plot_success_rates, summary_table
    from dockprepbench.analysis.stats import load_run

    runs = {_run_name(path): load_run(_results_path(path)) for path in args.runs}
    if len(runs) != len(args.runs):
        raise SystemExit("Run names must be unique; each run directory needs its own name")
    table = summary_table(runs)
    args.out.mkdir(parents=True, exist_ok=True)
    table.to_csv(args.out / "summary.csv", float_format="%.4f")
    written = [args.out / "summary.csv"]
    written += plot_success_rates(table, args.out / "success_rates")
    written += plot_rmsd_cdf(runs, args.out / "rmsd_cdf")
    print("\n".join(str(path) for path in written))
    return 0


def _run_name(path: Path) -> str:
    """
    Name a run after its directory (the config name), or the CSV's stem otherwise.

    Args:
        path (Path): Run directory or results CSV.

    Returns:
        str: Display name.
    """
    if path.is_dir():
        return path.name
    return path.parent.name if path.name == RESULTS_FILE else path.stem


def _results_path(path: Path) -> Path:
    """
    Resolve a run directory or results CSV to the CSV path.

    Args:
        path (Path): Run directory or CSV file.

    Returns:
        Path: Path to the results CSV.
    """
    return path / RESULTS_FILE if path.is_dir() else path


def _format_summary(path: Path) -> str:
    """
    Describe a run's success rates and failures.

    Args:
        path (Path): Results CSV.

    Returns:
        str: Multi-line summary.
    """
    from dockprepbench.analysis.stats import load_run, success_rate, successes

    run = load_run(path)
    failed = run[run["error"] != ""]
    lines = [
        f"{path}: {len(run)} complexes, {len(failed)} failed",
        f"  Top-1 success (RMSD < 2 A): {success_rate(successes(run))}",
        f"  PB-valid success:           {success_rate(successes(run, pb_valid=True))}",
    ]
    if len(failed):
        kinds = failed["error"].str.split(":").str[0].value_counts()
        lines.append("  failures: " + ", ".join(f"{k} {v}" for k, v in kinds.items()))
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    """
    Entry point for the `dockprepbench` command.

    Args:
        argv (list[str] | None): Arguments to parse; defaults to sys.argv.

    Returns:
        int: Process exit code.
    """
    args = build_parser().parse_args(argv)
    setup_logging(logging.DEBUG if args.verbose else logging.INFO)
    commands = {
        "run": run_command,
        "summarize": summarize_command,
        "compare": compare_command,
        "figures": figures_command,
        "report": report_command,
    }
    if args.command in commands:
        return commands[args.command](args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
