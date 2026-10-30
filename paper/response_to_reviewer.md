# Response to the reviewer

We thank the reviewer for a careful and constructive report. We agree with the central
point: the manuscript is strongest as a reproducible benchmark of preparation sensitivity in
redocking, and some claims went beyond that setting. We have narrowed those claims, added the
requested analyses and added an independent check whose analysis plan was committed before
any of its runs. Section and table numbers refer to the revised manuscript. Reviewer comments
are summarised in italics.

## Major comment 1: novelty and positioning

*Distinguish the study from earlier preparation-variability work (Corbeil et al., Sastry et
al.); explain what the paired design adds; present the agreement analysis as a separate
contribution; avoid implying that the importance of preparation is shown for the first
time.*

- **Introduction, new paragraph.** It states that the importance of preparation is not new,
  and that earlier work reported effects on aggregate success or enrichment, typically for
  one program. Our contribution is a common scale: preparation effects measured against
  engine differences and against seed-to-seed variation. The paragraph also explains why
  only a paired design shows which complexes a choice changes.
- **Contributions list.** The agreement analysis is now a separate, explicitly exploratory
  methodological item.

## Major comment 2: the crystal-ligand box

*The box is derived from the crystal ligand, so the 2 Å result is conditional on redocking
and does not show that 2 Å is optimal for prospective docking.*

- **Scope stated throughout.** "Redocking with the search space drawn from the crystal
  ligand" now appears in the abstract, introduction, discussion and conclusions.
- **Discussion, box recommendation reworded.** The tight-box result shows that a generously
  drawn box around a known site costs accuracy. It does *not* show that 2 Å is optimal for
  prospective docking, cross-docking or screening. We note that a tight box can cut off the
  binding mode, as our ligand-scaled box did (Section 3.9).
- **The independent Astex check (Section 3.11) supports the reviewer's caution.** The 2 Å gain
  recurred only for rDock (+17.6 points [8.2, 26.4]). For Vina, Smina and Gnina there was no
  detectable effect, and for Vina the interval (+1.2 [−2.9, 5.2]) excludes the benchmark-set
  gain. The geometric pattern did recur: the tight box removed wrong-site poses and larger
  boxes added them. But the Astex baselines rarely produced wrong-site poses for the Vina
  family, so there was little to fix. The abstract, discussion and conclusions now say the
  benefit depends on the engine and the dataset.
- **Cross-docking** remains future work and is named as such in the Limitations.

## Major comment 3: receptor preparation and cofactors

*Explain why heterogens are removed; quantify cofactor-containing complexes and their
contribution; clarify each engine's receptor representation; avoid attributing pipeline
differences to docking algorithms alone.*

- **Methods 2.3, rationale.** Most cofactors have no Meeko template. Keeping them for some
  engines only would confound the engine comparison. The removal is acknowledged as not
  neutral.
- **New Table 1: receptor representation per engine** (file format, conversion, atom typing).
- **"Engine" differences are now called complete-pipeline differences** throughout:
  Methods 2.3, Results 3.5 and Discussion.
- **New robustness analysis (Section 3.7).**
  - 244 of 428 complexes contain a cofactor or ion.
  - 7 to 24 correct top poses per engine fail PoseBusters only on a cofactor check, the most
    for Gnina.
  - On the 184 complexes without cofactors or ions, the engines keep their order and the 2 Å
    box raises success for all four: Vina +6.5 [1.2, 11.8], Smina +7.1 [2.1, 12.0], Gnina
    +2.2 [−2.8, 7.1], rDock +9.2 [3.4, 14.9] points.
  - The Limitations note that keeping relevant cofactors could change PB-valid rankings.

## Major comment 4: the agreement-based confidence signal

*It is a candidate signal; report n and prevalence per agreement level; assess added value
beyond the docking score; evaluate on held-out or independent data; the claims should match
the evidence.*

- **Wording.** It is now a "candidate confidence signal" in the section title, discussion
  and conclusions, explicitly not a validated or calibrated measure.
- **n and success rate at every agreement level** are in Figure 3, with the counts in
  `paper/tables/agreement_levels.csv`. Failed runs count as disagreement.
- **Held-out evaluation (new Table 8): 20 × 5-fold cross-validation over complexes.**
  - Held-out AUROC: agreement 0.81–0.88 against 0.64–0.72 for the docking score.
  - Lower Brier scores.
- **Added value beyond the score.**
  - Score and agreement together: held-out AUROC 0.86–0.91.
  - A likelihood-ratio test for adding agreement to the score gives p < 10⁻³⁰ for every
    engine.
- **Missing poses** affect only 16 Vina and 2 rDock complexes. Excluding them changes Vina's
  AUROC from 0.857 to 0.847.
- **Independent check on the Astex Diverse Set (new Section 3.11).**
  - 85 complexes, no PDB entry in common with the benchmark set.
  - Four hypotheses, their tests and decision rules were committed to the repository
    (`docs/astex_validation_plan.md`) before any run.
  - **Agreement replicated for three of four engines.** It beat the docking score for Vina
    (AUROC 0.82 vs 0.66; difference 0.16 [0.04, 0.29]), Smina (0.80 vs 0.62) and rDock
    (0.93 vs 0.78). For Gnina the difference was similar in size but its interval included
    zero (0.14 [−0.01, 0.28]; only 17 incorrect poses), and agreement did not beat Gnina's
    CNN score.
  - **The pair chosen on the benchmark set** also beat the score out of sample for Vina,
    Smina and rDock, not for Gnina.
  - **Calibration held for all four engines.** When all five arms agreed, the pose was correct
    for 87–91% of complexes, against 81–93% on the benchmark set.
  - We report these results whatever their direction, as the plan required, including the
    box result below, which did not replicate.
- **Shared systematic errors.** We agree that the arms share receptor, scoring and much of
  the search, so agreement can reflect shared errors. The seed control (0.67) shows the
  signal is not mere run-to-run noise. Shared errors are why we describe the signal as
  candidate and test it on an independent set.

## Major comment 5: statistics

- **5.1 Paired intervals.** Every arm-vs-baseline difference now has a 95% interval for paired
  proportions (Newcombe's hybrid score interval, checked against a paired bootstrap). The
  intervals are in the text, Table 3 and SI Tables S1–S4.
- **5.2 Hypothesis families.** A new paragraph in Methods 2.7 separates the confirmatory
  analyses from the exploratory ones.
  - *Confirmatory:* the five primary arms per engine and metric, fixed before the runs.
  - *Exploratory:* the follow-up arms, ligand properties, agreement and robustness analyses.
  - It also explains why the overlapping follow-up arms are corrected separately.
- **5.3 Engine comparisons** are described as differences between complete pipelines (see
  comment 3).
- **5.4 Seed comparison.** The "up to 4.9 times" ratios are now explicitly descriptive in
  Results 3.3.

## Major comment 6: protonation and conformer

- **Protonation** is now described as the supplied state against the single dominant state
  predicted at pH 7.4. We state that protomer and tautomer ensembles and receptor protonation
  were not tested (Methods 2.4; Limitations).
- **Conformer** is now described as two starting-structure protocols that differ in
  minimisation force field (UFF vs MMFF94) as well as coordinates. Its effect cannot be
  assigned to the coordinates alone (Methods 2.4; Discussion).

## Major comment 7: failure-mode classes

- **Classes as descriptions, not diagnoses.** They are now defined as geometric descriptions,
  not diagnoses of why docking failed, with examples of where they mislead (Methods 2.6).
  Results and Discussion use "describe" and "suggest" instead of "explain".
- **Threshold sensitivity (new Table 5).** Wrong-site counts at 2–6 Å: at every threshold
  and for every engine, the 2 Å box gives fewer wrong-site poses than the baseline and the
  8 Å box more. The box-size conclusion does not depend on the 4 Å choice.

## Major comment 8: generalisability

- The 308-complex and cofactor analyses are now labelled within-dataset robustness checks
  (Results 3.6, 3.7; Limitations).
- The Astex check is described as one small independent test, not a general validation.
- Cross-docking, apo and predicted structures, and independently defined sites are named as
  future tests.

## Minor comments

1. **Range vs paired effect.** Results 3.5 now says the spread across protocols is a
   best-to-worst range, not the effect of one choice.
2. **33.2 points.** Now described as a difference between complete tested pipelines, not
   between neural-network and classical methods.
3. **Top-1 vs PB-valid.** Each is defined in Methods 2.6, which now also states that
   "success" without qualification always means top-1 success. We checked every use in
   Results, Discussion, Abstract and Conclusions: PB-valid success is always named as such.
4. **Missing poses in the agreement analysis.** Their handling is stated in Methods and the
   analysis is repeated without them (see comment 4).
5. **Selection optimism of the best pair.** Choosing the pair inside each training fold gives
   held-out AUROC 0.79–0.87, about 0.02 below the in-sample values. The pre-selected pair is
   tested out of sample on Astex (H3).
6. **Reproducibility of the numbers.** 312 of the 369 decimal numbers in Results are now generated directly from the analysis output tables (`\val{}` references, regenerated with the tables), so the text cannot drift from the data. The rest are settings (pH, versions), simple derived ranges and ratios, or a few one-off statistics. Every follow-up analysis that was previously computed ad hoc (protonation x box, ligand-sized box, engine comparisons, sampling counts, comparison with the published PoseBusters run) now has its own module and table.
7. **Zenodo DOI.** \[To be added on release.\]
8. **Authors and acknowledgements.** \[To be completed.\]
9. **Paired nature in the figures.** New SI Figure (gained and lost complexes per arm and
   engine) shows each net change as the difference of two often much larger counts, with the
   seed-to-seed level marked for Vina. The explorer also shows, for every arm, which
   complexes were gained and lost.
