# Independent validation on the Astex Diverse Set: analysis plan

Written and committed before any Astex docking run, so that the hypotheses, analyses and
decision rules below are fixed in advance. The commit timestamp is the record. Anything
analysed beyond this plan will be reported as exploratory.

## Why

The PoseBusters analysis found two results that a reviewer asked to see confirmed on
independent data:

1. **Box padding.** A 2 Å box raised top-1 success for all four engines, significantly for
   Vina, Smina and rDock.
2. **Agreement across preparations as a confidence signal.** The number of preparation arms
   whose top-1 pose lies within 2 Å of the baseline pose identified correct baseline poses
   better than the docking score (AUROC 0.83–0.90 against 0.64–0.72). The analysis was
   exploratory. The pair "2 Å box + RDKit conformer" was the best two-arm subset for every
   engine, but it was selected on the same data.

## Data

The Astex Diverse Set (85 complexes) from the same PoseBusters Zenodo release
(record 8278563). It shares no PDB entry with the 428-complex benchmark set. Every complex
is included; preparation or docking failures count as failures, and a missing pose counts
as disagreement, as in the main study. No complex is excluded after the runs.

## Runs

The six protocols of the main study, unchanged: baseline (4 Å padding, supplied ligand),
box padding 2, 6 and 8 Å, protonation at pH 7.4 (Dimorphite-DL) and the RDKit start
conformer. Each runs with Vina, Smina, Gnina and rDock, using the same engine settings and
seed as the main study (`dockprepbench run --dataset astex`). Results go to
`results/astex/<engine>/<arm>/`.

## Hypotheses and tests

**H1 (box).** The 2 Å box raises top-1 success over the baseline.
Test: for each engine, the paired difference with a 95% interval, discordant counts, and
the exact McNemar test. Holm correction over the four engines.
H1 is supported for an engine if the difference is positive and Holm p < 0.05.
With 85 complexes the power is low, so the direction and the interval are reported for
every engine whatever the p-value.

**H2 (agreement, five arms).** Agreement with the five arms (0–5) has a higher AUROC for
baseline top-1 success than the docking score (negated). For Gnina it is also compared with
the CNN score.
Test: the paired difference in AUROC, with a 95% interval from 2,000 bootstrap resamples of
complexes (seed 0).
H2 is supported for an engine if the interval excludes zero.

**H3 (agreement, the pre-selected pair).** Agreement with the two arms chosen on PoseBusters,
the 2 Å box and the RDKit conformer (0–2), has a higher AUROC than the docking score.
Same test and rule as H2. This is the out-of-sample test of the pair selection.

**H4 (calibration).** When all five arms agree, the baseline pose is correct about as often
as on PoseBusters.
Test: the success rate among unanimous complexes with a Wilson interval, compared with the
PoseBusters rate (Vina 83.7%, Smina 81.3%, Gnina 92.8%, rDock 91.1%).
This is descriptive: overlapping intervals are read as consistent.

## Expected caveats, stated in advance

- **High baseline success.** Astex is an older, easier set, so baseline success is likely
  high. Few incorrect poses leave little room for a confidence signal to discriminate, and
  the AUROC intervals will be wide. If an engine has fewer than 10 incorrect baseline
  poses, its AUROC will be reported but not tested.
- **One small set.** 85 complexes is one independent check, not a general validation.
  Cross-docking remains untested.

## Reporting

All four hypotheses are reported for all four engines, whatever the outcome, in a
validation section of the paper and in the response to reviewers.
