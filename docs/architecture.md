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
(the `dockprepbench` command).

## Repository layout

```
src/dockprepbench/   library code
tests/               pytest suite, mirroring the package
docs/                documentation
scripts/             one-off utilities (not imported by the package)
data/                downloaded inputs (git-ignored)
results/             generated outputs (git-ignored)
```
