"""Command-line interface: `dockprepbench`."""

import argparse
import logging

from dockprepbench import __version__
from dockprepbench.log import setup_logging


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="dockprepbench",
        description="Benchmark how preparation choices affect classical docking.",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable debug logging.")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    setup_logging(logging.DEBUG if args.verbose else logging.INFO)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
