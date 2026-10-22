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
