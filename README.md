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

On a Linux machine or shared server, one script installs everything in user space (uv, the
Python environment, Smina, Gnina with cuDNN 9, rDock and Open Babel), keeps caches out of your
home directory and finishes by running the tests:

```bash
scripts/setup.sh --tools-dir /path/with/space/tools   # add --no-gnina to skip the 2 GB Gnina binary
source /path/with/space/tools/env.sh
uv run dockprepbench run --engine smina --workers 16 --config configs/*.toml
```

Engines: `--engine vina|smina|gnina|rdock`; Gnina spreads workers over GPUs with `--gpus 0,1,2`.

## Usage

Docking runs on Linux (see Installation). Redock the PoseBusters benchmark set with Vina
(downloads the data to `data/raw` on first run):

```bash
uv run dockprepbench run --workers 8                                   # baseline protocol
uv run dockprepbench run --config configs/box_padding_2.toml --workers 8
```

Each preparation protocol is a TOML file in [`configs/`](configs/) that states only what it
changes from the baseline. A run writes `results.csv` (one row per complex: top-1 RMSD,
PoseBusters validity, score, runtime, any error) and `run.json` (config, engine settings,
version) to `results/<engine>/<config name>/`. Prepared receptors are cached in
`data/prepared/`, shared between configs that prepare the receptor identically.

```bash
uv run dockprepbench summarize results/vina/baseline
uv run dockprepbench compare results/vina/baseline results/vina/box_padding_2
```

`summarize` reports top-1 and PB-valid success rates with 95% Wilson intervals. `compare`
pairs two runs on their shared complexes, reports an exact McNemar test and lists the
complexes whose outcome changed. Use `--limit N` on `run` for a quick check.

## Project layout

```
src/dockprepbench/
├── datasets/      benchmark download, curation and splits
├── preparation/   receptor and ligand preparation variants
├── engines/       docking engine wrappers
├── evaluation/    RMSD, PoseBusters validity, enrichment
├── analysis/      aggregation, statistics, figures
├── benchmark.py   run an engine over a dataset and score poses
├── io.py          structure file I/O
├── log.py         logging setup
└── cli.py         `dockprepbench` command
configs/           preparation protocols (TOML), one per experiment arm
tests/             pytest suite mirroring the package
docs/              documentation
scripts/           one-off utilities
```

See [docs/methods.md](docs/methods.md) for the benchmark protocol, [docs/architecture.md](docs/architecture.md)
for how the stages fit together, and
[CONTRIBUTING.md](CONTRIBUTING.md) for the development workflow.

## Roadmap

- [x] Project foundation: packaging, tests, CI
- [ ] Vina end-to-end on a small complex set with RMSD scoring
- [ ] All engines + PoseBusters validity checks on the full set
- [ ] Preparation-variant experiments
- [ ] Figures, write-up and preprint

## License

MIT
