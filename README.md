# DockPrepBench

[![CI](https://github.com/David-Alobo/DockPrepBench/actions/workflows/ci.yml/badge.svg)](https://github.com/David-Alobo/DockPrepBench/actions/workflows/ci.yml)

**How much does preparation matter in classical docking?**

DockPrepBench measures how much protein and ligand preparation choices change docking outcomes,
compared with the choice of docking engine.

Docking benchmarks usually compare engines under a single, fixed preparation protocol. In practice,
decisions such as protonation state, tautomer enumeration, retained waters, holo vs. apo receptor and
binding-box definition can shift results as much as switching tools. This project quantifies that
effect in a reproducible, hardware-agnostic pipeline.

> Status: early development. Results will be added as milestones are completed.

## Scope

| | |
|---|---|
| **Engines** | AutoDock Vina, Smina, Gnina, rDock/CmDock |
| **Datasets** | PoseBusters benchmark set; similarity-split PLINDER subset; virtual screening subset (planned) |
| **Metrics** | Top-1 success (RMSD < 2 Å), PoseBusters physical validity, enrichment (EF, BEDROC) |

## Installation

Requires [uv](https://docs.astral.sh/uv/).

```bash
uv sync
uv run pytest
```

The pipeline runs on CPU. A GPU is optional and only speeds up Gnina.

Docking engines run on **Linux x86_64** only (the `vina` Python bindings ship Linux wheels only, and
Gnina and rDock are Linux-only). On Windows or macOS the library installs and the unit tests run, while
engine tests are skipped; use WSL, a Linux machine or CI to dock.

## Project layout

```
src/dockprepbench/
├── datasets/      benchmark download, curation and splits
├── preparation/   receptor and ligand preparation variants
├── engines/       docking engine wrappers
├── evaluation/    RMSD, PoseBusters validity, enrichment
├── analysis/      aggregation, statistics, figures
├── io.py          structure file I/O
├── log.py         logging setup
└── cli.py         `dockprepbench` command
tests/             pytest suite mirroring the package
docs/              documentation
scripts/           one-off utilities
```

See [docs/architecture.md](docs/architecture.md) for how the stages fit together and
[CONTRIBUTING.md](CONTRIBUTING.md) for the development workflow.

## Roadmap

- [x] Project foundation: packaging, tests, CI
- [ ] Vina end-to-end on a small complex set with RMSD scoring
- [ ] All engines + PoseBusters validity checks on the full set
- [ ] Preparation-variant experiments
- [ ] Figures, write-up and preprint

## License

MIT
