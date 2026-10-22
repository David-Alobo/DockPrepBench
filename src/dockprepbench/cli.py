"""Command-line interface: `dockprepbench`."""

import argparse
import json
import logging
from datetime import UTC, datetime
from pathlib import Path

from dockprepbench import __version__
from dockprepbench.log import setup_logging

logger = logging.getLogger(__name__)

RESULTS_FILE = "results.csv"


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
    run.add_argument(
        "--config", type=Path, help="Preparation config (TOML). Defaults to the baseline."
    )
    run.add_argument("--data-dir", type=Path, default=Path("data/raw"))
    run.add_argument("--work-dir", type=Path, default=Path("data/prepared"))
    run.add_argument(
        "--out", type=Path, help="Run directory. Defaults to results/<engine>/<config name>."
    )
    run.add_argument("--limit", type=int, help="Dock only the first N complexes.")
    run.add_argument("--workers", type=int, default=1, help="Complexes docked in parallel.")
    run.add_argument("--exhaustiveness", type=int, default=8)

    summarize = commands.add_parser("summarize", help="Success rates of one run.")
    summarize.add_argument("run", type=Path, help="Run directory or results CSV.")

    compare = commands.add_parser("compare", help="Paired comparison of two runs.")
    compare.add_argument("run_a", type=Path, help="Run directory or results CSV.")
    compare.add_argument("run_b", type=Path, help="Run directory or results CSV.")
    compare.add_argument(
        "--pb-valid", action="store_true", help="Require PB-valid poses for success."
    )
    return parser


def run_command(args: argparse.Namespace) -> int:
    """
    Execute the `run` subcommand.

    Writes `results.csv` and `run.json` (config, engine settings, version, time) to the
    run directory.

    Args:
        args (argparse.Namespace): Parsed `run` arguments.

    Returns:
        int: Process exit code.
    """
    from dockprepbench.benchmark import run_benchmark
    from dockprepbench.datasets.posebusters import download, load_complexes
    from dockprepbench.engines.vina import VinaEngine
    from dockprepbench.preparation.config import PrepConfig, load_config

    config = load_config(args.config) if args.config else PrepConfig()
    out = args.out or Path("results") / args.engine / config.name
    complexes = load_complexes(download(args.data_dir))[: args.limit]
    # One Vina thread per complex when complexes already run in parallel.
    engine = VinaEngine(exhaustiveness=args.exhaustiveness, cpu=1 if args.workers > 1 else 0)
    logger.info("Docking %d complexes with %s, config %s", len(complexes), engine.name, config)

    started = datetime.now(UTC)
    results = run_benchmark(complexes, engine, config, args.work_dir, workers=args.workers)
    out.mkdir(parents=True, exist_ok=True)
    results.to_csv(out / RESULTS_FILE, index=False)
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

    logger.info("Results written to %s", out)
    print(_format_summary(out / RESULTS_FILE))
    return 0


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
    Execute the `compare` subcommand.

    Args:
        args (argparse.Namespace): Parsed `compare` arguments.

    Returns:
        int: Process exit code.
    """
    from dockprepbench.analysis.stats import compare, load_run, success_rate, successes

    path_a, path_b = _results_path(args.run_a), _results_path(args.run_b)
    hits_a = successes(load_run(path_a), pb_valid=args.pb_valid)
    hits_b = successes(load_run(path_b), pb_valid=args.pb_valid)
    result = compare(hits_a, hits_b)
    shared = hits_a.index.intersection(hits_b.index)

    metric = "PB-valid success" if args.pb_valid else "Top-1 success"
    lines = [
        f"{metric} on {result.n} shared complexes",
        f"  A {path_a}: {success_rate(hits_a[shared])}",
        f"  B {path_b}: {success_rate(hits_b[shared])}",
        f"  both {result.both}, neither {result.neither}, "
        f"only A {len(result.only_a)}, only B {len(result.only_b)}",
        f"  McNemar exact p = {result.p_value:.3g}",
    ]
    if result.only_a:
        lines.append(f"  only A: {' '.join(result.only_a)}")
    if result.only_b:
        lines.append(f"  only B: {' '.join(result.only_b)}")
    print("\n".join(lines))
    return 0


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
    commands = {"run": run_command, "summarize": summarize_command, "compare": compare_command}
    if args.command in commands:
        return commands[args.command](args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
