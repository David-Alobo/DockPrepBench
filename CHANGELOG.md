# Changelog

All notable changes to this project are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the project uses
[Semantic Versioning](https://semver.org/).

## [Unreleased]

### Added
- Package layout (`dockprepbench`) with `datasets`, `preparation`, `engines`, `evaluation` and `analysis` subpackages.
- `dockprepbench` command-line entry point.
- uv-managed environment, pytest suite, ruff linting, pre-commit hooks and cross-platform CI.
- `evaluation.rmsd`: symmetry-corrected, in-place ligand RMSD and top-1 success rate (< 2 Å).
- `io.load_sdf_file` for reading ranked docked poses.
- `engines.base`: `DockingEngine` interface with `Box` and `DockingResult`.
- `engines.vina.VinaEngine`: AutoDock Vina via its Python bindings (Linux x86_64 only), with a fixed default seed.
- `engines.pdbqt`: meeko-based receptor/ligand PDBQT conversion and pose read-back with bond orders preserved.
- `preparation.box.box_from_ligand`: search box around a reference ligand with 4 Å padding.
