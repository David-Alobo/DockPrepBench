# Changelog

All notable changes to this project are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the project uses
[Semantic Versioning](https://semver.org/).

## [Unreleased]

### Added
- Every Results number traceable to data: new tables for engine comparisons, sampling counts, the comparison with the published PoseBusters Vina run, the follow-up arms (`analysis.followup`) and pose-change rates; 188 more numbers migrated to `\val{}` (312 of 369), compiled text unchanged. SI figure of gained and lost complexes per arm (`plot_flips`); Methods defines unqualified success as top-1.
- Results numbers taken from the data: 124 numbers in `results.tex` are `\val{}` references to values generated from the result tables (`paper/tables/numbers.tex`, spec `paper/numbers.toml`, generated with provenum); the compiled text is unchanged.
- Independent check on the Astex Diverse Set (85 complexes; plan in `docs/astex_validation_plan.md`, `analysis.astex`): agreement beat the docking score for Vina, Smina and rDock with calibration held for all four engines; the 2 Å box gain recurred only for rDock. Paper section 3.11, Table 10; abstract, discussion, conclusions and response letter updated.
- `analysis.agreement`: pose agreement across preparation arms, engines and seeds as a confidence signal (AUROC with paired bootstrap intervals), the consensus pose, and the per-complex oracle; `analysis.figures.plot_agreement_levels`; how many arms the signal needs (`subset_table`); paper section 3.9 and Figure 3.
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
- `preparation.box.box_from_gyration` and the `box_size` config field: a cubic box 2.857 times
  the docked ligand's radius of gyration (Feinstein and Brylinski 2015); arm `box_gyration`.
- Interaction arms `ligand_dimorphite_box_2`, `_6` and `_8`: ligand protonation combined with
  each box padding.

### Changed
- Pose RMSD matches atoms on the charge-free heavy-atom skeleton, so poses docked in a different
  protonation state from the crystal ligand still score.
- Receptor PDBQT conversion drops residues that fail template matching when they lie more
  than 8 Å outside the search box, and still fails for residues near the box.
