# Skills demonstrated

DockPrepBench is a research-software project that measures how much protein and ligand
preparation choices change molecular docking outcomes, compared with the choice of
docking engine. This page summarises the skills it demonstrates, each with a pointer to
where it can be seen in the repository.

## Computational chemistry

- **Protein structure preparation.** Repairing crystal structures (missing atoms,
  nonstandard residues, chain fragments), assigning protonation states at a given pH, and
  relaxing added atoms while keeping experimental coordinates fixed.
  → [`preparation/receptor.py`](../src/dockprepbench/preparation/receptor.py)
- **Ligand preparation.** pH-dependent protonation that preserves 3D geometry, and
  seeded conformer generation (ETKDGv3, MMFF94).
  → [`preparation/ligand.py`](../src/dockprepbench/preparation/ligand.py)
- **Molecular docking.** AutoDock Vina through its Python API: search-space definition,
  PDBQT conversion with meeko, and pose recovery with bond orders intact.
  → [`engines/`](../src/dockprepbench/engines/)
- **Pose evaluation.** Symmetry-corrected RMSD without superposition, robust to
  protonation-state differences, plus physical-validity checks with PoseBusters.
  → [`evaluation/`](../src/dockprepbench/evaluation/)
- **Cheminformatics with RDKit.** Sanitisation, substructure matching, stereochemistry
  from 3D coordinates, and charge-agnostic atom mapping.
  → [`chem.py`](../src/dockprepbench/chem.py)

## Research methods and statistics

- **Experimental design.** A one-factor-at-a-time benchmark in which each experiment arm
  changes a single preparation choice from a baseline, with every arm defined in a
  versioned config file. → [`configs/`](../configs/), [methods](methods.md)
- **Paired statistics.** Comparing arms on the same complexes with exact McNemar
  tests, and reporting success rates with Wilson score intervals.
  → [`analysis/stats.py`](../src/dockprepbench/analysis/stats.py)
- **Benchmark integrity.** Checksum-verified data, explicit exclusion criteria with a
  recorded reason for every failure, and full-dataset validation sweeps before
  expensive runs. These sweeps found preparation failures that unit tests missed.
- **Debugging scientific software.** Traced intermittent failures to their root causes:
  unseeded random placement in OpenMM, order-dependent multithreaded force sums, and
  distance-based bond perception. Each fix was confirmed by measurement, for example
  making receptor preparation byte-for-byte reproducible.

## Research software engineering

- **Python packaging and tooling.** An installable package with a command-line
  interface, managed with uv and checked with ruff and pre-commit.
  → [`pyproject.toml`](../pyproject.toml), [`cli.py`](../src/dockprepbench/cli.py)
- **Testing.** A pytest suite of more than 80 tests, including regression tests for
  every bug found and stand-in implementations for platform-specific components.
  → [`tests/`](../tests/)
- **Continuous integration.** GitHub Actions on Linux, Windows and macOS, with test
  failures shown inline on pull requests. → [`.github/workflows/`](../.github/workflows/)
- **Version control.** Feature branches, focused pull requests and descriptive commit
  history.
- **Reproducibility.** Pinned dependencies, deterministic preparation and docking, and
  a record of the configuration, engine settings and version written next to every
  result.
- **Cross-platform work.** Developed on Windows; docking runs on Linux under WSL
  because some dependencies are Linux-only. Parallel processing across CPU cores.

## Scientific communication

- **Publication figures.** Colour-blind-safe palettes validated by script, consistent
  colour across figures, vector output, and a study-design schematic generated from code.
  → [`analysis/figures.py`](../src/dockprepbench/analysis/figures.py),
  [`figures/concept.svg`](figures/concept.svg)
- **Technical writing.** A methods section kept in step with the code, stating known
  limitations. → [methods](methods.md), [architecture](architecture.md)

## In progress

Planned stages of the project, not yet demonstrated:

- Results across all experiment arms and their interpretation
- Additional docking engines (Smina, Gnina, rDock) and a comparison of preparation
  and engine effects
- Running at scale on a high-performance computing cluster
- Preprint write-up
