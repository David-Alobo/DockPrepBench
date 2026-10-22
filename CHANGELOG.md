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
- `datasets.posebusters`: checksum-verified download and loading of the PoseBusters benchmark set.
- `preparation.receptor.protonate_receptor`: PDBFixer repair and protonation at pH 7.4, heterogens removed;
  added atoms are relaxed with Amber14 while crystal atoms stay fixed.
- `benchmark`: redock a set of complexes in parallel and record top-1 RMSD, score and failures per complex.
- `dockprepbench run` command for a PoseBusters redocking run.
- `evaluation.validity`: PoseBusters physical-validity checks on the top pose; runs report
  PB-valid success (RMSD < 2 Å and every check passed) alongside plain top-1 success.
- `preparation.config`: named preparation protocols (`PrepConfig`) loaded from TOML files in
  `configs/` (baseline and box-padding 2/6/8 Å); `dockprepbench run --config` writes results and
  `run.json` to `results/<engine>/<config>/`, caching receptors per receptor-relevant setting.
- `analysis.stats`: success rates with Wilson intervals and paired McNemar comparison of runs;
  `dockprepbench summarize` and `dockprepbench compare` commands.
- `preparation.ligand`: ligand protonation at the config pH (Dimorphite-DL, keeping heavy-atom
  coordinates) and RDKit conformer regeneration (ETKDGv3 + MMFF, seeded); config fields
  `ligand_protonation` and `ligand_conformer`, with `ligand_dimorphite` and
  `ligand_rdkit_conformer` configs.
- `dockprepbench run --config a.toml b.toml ...` runs several configs in turn, skipping finished
  ones (`--force` reruns); `dockprepbench figures` writes a success-rate forest plot, RMSD curves
  and a summary table.
- `docs/methods.md`: the benchmark protocol as implemented.
- Results record `best_rmsd` and `best_rank` (best pose among all returned, to separate scoring
  from sampling failures) and `top1_centroid_shift` (to separate wrong-site from wrong-mode poses).
- `dockprepbench compare BASELINE ARM [ARM ...]` compares several arms at once and reports
  Holm-adjusted p-values (`analysis.stats.holm_adjust`).
- Runs save every returned pose to `results/<engine>/<arm>/poses/<id>.sdf`, in rank order with
  `rank`, `score` and `rmsd` properties, for inspection and the planned results explorer.

### Changed
- Pose RMSD matches atoms on the charge-free heavy-atom skeleton, so poses docked in a different
  protonation state from the crystal ligand still score.
- Receptor PDBQT conversion drops residues that fail template matching when they lie more
  than 8 Å outside the search box, and still fails for residues near the box.
