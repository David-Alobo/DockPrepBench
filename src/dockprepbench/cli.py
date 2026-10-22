"""Command-line interface: `dockprepbench`."""

import argparse
import logging
from pathlib import Path

from dockprepbench import __version__
from dockprepbench.log import setup_logging

logger = logging.getLogger(__name__)


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
    run.add_argument("--engine", choices=["vina"], default="vina")
    run.add_argument("--data-dir", type=Path, default=Path("data/raw"))
    run.add_argument("--work-dir", type=Path, default=Path("data/prepared"))
    run.add_argument("--out", type=Path, required=True, help="Output CSV, one row per complex.")
    run.add_argument("--limit", type=int, help="Dock only the first N complexes.")
    run.add_argument("--workers", type=int, default=1, help="Complexes docked in parallel.")
    run.add_argument("--exhaustiveness", type=int, default=8)
    return parser


def run_command(args: argparse.Namespace) -> int:
    """
    Execute the `run` subcommand.

    Args:
        args (argparse.Namespace): Parsed `run` arguments.

    Returns:
        int: Process exit code.
    """
    from dockprepbench.benchmark import run_benchmark
    from dockprepbench.datasets.posebusters import download, load_complexes
    from dockprepbench.engines.vina import VinaEngine
    from dockprepbench.evaluation.rmsd import top1_success_rate

    complexes = load_complexes(download(args.data_dir))[: args.limit]
    # One Vina thread per complex when complexes already run in parallel.
    engine = VinaEngine(exhaustiveness=args.exhaustiveness, cpu=1 if args.workers > 1 else 0)
    logger.info("Docking %d complexes with %s", len(complexes), engine.name)

    results = run_benchmark(complexes, engine, args.work_dir, workers=args.workers)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    results.to_csv(args.out, index=False)

    failed = (results["error"] != "").sum()
    logger.info(
        "Top-1 success: %.1f%% (%d complexes, %d failed). Results: %s",
        100 * top1_success_rate(results["top1_rmsd"]),
        len(results),
        failed,
        args.out,
    )
    return 0


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
    if args.command == "run":
        return run_command(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
