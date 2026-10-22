/**
 * A small explorer dataset whose expected numbers can be worked out by hand.
 *
 * 25 complexes, C00_LIG ... C24_LIG. Baseline: C00–C09 valid successes, C10–C14 scoring
 * misses, C15–C19 wrong site, C20–C24 wrong binding mode (10/25 = 40% success).
 * box_padding_2: like the baseline, except C00 is lost (wrong mode) and C10–C14 are gained
 * as valid successes (14/25 = 56%; 5 gained, 1 lost). C05 is a correct pose that fails
 * PoseBusters in both arms. ligand_dimorphite equals the baseline.
 * Every fifth complex (C00, C05, ...) has a cofactor.
 */

const ok = (valid = true) => ({ rmsd: 1.0, best: 1.0, valid, shift: 0.5 });
const scoringMiss = () => ({ rmsd: 4.0, best: 1.5, valid: true, shift: 2.0 });
const wrongSite = () => ({ rmsd: 9.0, best: 6.0, valid: true, shift: 7.0 });
const wrongMode = () => ({ rmsd: 5.0, best: 4.0, valid: true, shift: 1.5 });

function baselineOutcome(i) {
  if (i < 10) return ok(i !== 5);
  if (i < 15) return scoringMiss();
  if (i < 20) return wrongSite();
  return wrongMode();
}

function boxOutcome(i) {
  if (i === 0) return wrongMode();
  if (i >= 10 && i < 15) return ok();
  return baselineOutcome(i);
}

export function makeData() {
  const complexes = Array.from({ length: 25 }, (_, i) => ({
    id: `C${String(i).padStart(2, "0")}_LIG`,
    heavy: 10 + i,
    rotb: i % 7,
    mw: 150 + 10 * i,
    cofactors: i % 5 === 0 ? ["HEM"] : [],
    svg: "<svg></svg>",
    arms: {
      baseline: baselineOutcome(i),
      box_padding_2: boxOutcome(i),
      ligand_dimorphite: baselineOutcome(i),
    },
  }));
  return {
    demo: true,
    shareUrl: "https://example.org/explorer",
    dataset: "Fixture set",
    engine: "Test engine",
    arms: [
      { key: "baseline", label: "Baseline", detail: "reference" },
      { key: "box_padding_2", label: "Box 2 Å", detail: "Box padding 2 Å" },
      { key: "ligand_dimorphite", label: "Protonated ligand", detail: "Ligand protonated" },
    ],
    complexes,
    real: {},
  };
}
