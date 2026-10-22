# Architecture

DockPrepBench is a pipeline. Each stage is a subpackage of `dockprepbench`:

```
datasets ──► preparation ──► engines ──► evaluation ──► analysis
```

| Subpackage | Responsibility |
|---|---|
| `datasets` | Download and curate benchmark sets (PoseBusters, PLINDER) and define splits. |
| `preparation` | Generate receptor and ligand variants: protonation, tautomers, retained waters, holo/apo, box definition. |
| `engines` | Run Vina, Smina, Gnina and rDock through a common interface. |
| `evaluation` | Score poses (RMSD, PoseBusters validity) and screens (EF, BEDROC). |
| `analysis` | Aggregate results, run statistics and make figures. |

Shared helpers live at the package root: `io` (structure files), `log` (logging) and `cli`
(the `dockprepbench` command). `benchmark` ties the stages together for one run: it prepares
each receptor (cached under `data/prepared`), docks the ligand from its generated start
conformer into a box around the crystal ligand, and scores the top-ranked pose.

## Repository layout

```
src/dockprepbench/   library code
tests/               pytest suite, mirroring the package
docs/                documentation
scripts/             one-off utilities (not imported by the package)
data/                downloaded inputs (git-ignored)
results/             generated outputs (git-ignored)
```
