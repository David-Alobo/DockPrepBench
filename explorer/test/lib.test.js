import { createRequire } from "node:module";
import { describe, expect, it } from "vitest";

import { makeData } from "./fixture.js";

const require = createRequire(import.meta.url);
const L = require("../lib.js");

const outcome = (o) => ({ rmsd: 5, best: 5, valid: true, shift: 0, ...o });

describe("classify", () => {
  it.each([
    [{ rmsd: 1.2, valid: true }, "valid"],
    [{ rmsd: 1.2, valid: false }, "invalid"],
    [{ rmsd: 3.0, best: 1.9 }, "scoring"],
    [{ rmsd: 3.0, best: 2.0, shift: 4.1 }, "site"],
    [{ rmsd: 3.0, best: 2.0, shift: 4.0 }, "mode"],
  ])("%o is %s", (o, expected) => {
    expect(L.classify(outcome(o))).toBe(expected);
  });

  it("treats 2 Å as a miss, as the Python evaluation does", () => {
    expect(L.classify(outcome({ rmsd: 2.0, best: 2.0 }))).not.toMatch(/valid/);
  });

  it("classifies a complex that raised an error as failed, before any RMSD rule", () => {
    expect(L.classify({ error: "RuntimeError: x" })).toBe("failed");
    expect(L.classify(outcome({ rmsd: 1, error: "x" }))).toBe("failed");
  });

  it("accepts a different site threshold", () => {
    expect(L.classify(outcome({ rmsd: 3, best: 3, shift: 3.5 }), 3)).toBe("site");
  });
});

describe("succeeds", () => {
  it("top-1 ignores PoseBusters validity; valid requires it", () => {
    const correctButInvalid = outcome({ rmsd: 1, valid: false });
    expect(L.succeeds(correctButInvalid, "top1")).toBe(true);
    expect(L.succeeds(correctButInvalid, "valid")).toBe(false);
    expect(L.succeeds(outcome({ rmsd: 2.5 }), "top1")).toBe(false);
  });

  it("never counts a failed complex as a success", () => {
    // A failed complex has no RMSD; `null < 2` is true in JavaScript, so this must be explicit.
    expect(L.succeeds({ error: "x", rmsd: null }, "top1")).toBe(false);
    expect(L.succeeds({ error: "x" }, "valid")).toBe(false);
  });
});

describe("wilson", () => {
  it("matches the reference interval used by the Python tests (8/10)", () => {
    const [lo, hi] = L.wilson(8, 10);
    expect(lo).toBeCloseTo(0.49, 3);
    expect(hi).toBeCloseTo(0.943, 3);
  });

  it("matches scipy-backed Python output for 163/428", () => {
    const [lo, hi] = L.wilson(163, 428);
    expect(lo).toBeCloseTo(0.33608953276238807, 12);
    expect(hi).toBeCloseTo(0.42771267307827054, 12);
  });

  it.each([[0, 5], [5, 5]])("stays inside [0, 1] for %i/%i", (k, n) => {
    const [lo, hi] = L.wilson(k, n);
    expect(lo).toBeGreaterThanOrEqual(0);
    expect(hi).toBeLessThanOrEqual(1);
    expect(lo).toBeLessThanOrEqual(k / n);
    expect(hi).toBeGreaterThanOrEqual(k / n);
  });

  it("rejects impossible counts", () => {
    expect(() => L.wilson(3, 0)).toThrow(RangeError);
    expect(() => L.wilson(6, 5)).toThrow(RangeError);
  });
});

describe("mcnemar", () => {
  it("equals the exact binomial tail for one-sided discordance", () => {
    expect(L.mcnemar(12, 0)).toBeCloseTo(2 * 0.5 ** 12, 15);
  });

  it("matches scipy.stats.binomtest for 75 vs 19", () => {
    expect(L.mcnemar(75, 19)).toBeCloseTo(4.8181440754468005e-9, 20);
  });

  it("is 1 with no or balanced discordance, and symmetric", () => {
    expect(L.mcnemar(0, 0)).toBe(1);
    expect(L.mcnemar(3, 3)).toBe(1);
    expect(L.mcnemar(4, 9)).toBeCloseTo(L.mcnemar(9, 4), 15);
  });
});

describe("rates and comparisons on the fixture", () => {
  const { complexes, arms } = makeData();

  it("computes success rates for both metrics", () => {
    expect(L.rate(complexes, "baseline", "top1")).toMatchObject({ k: 10, n: 25, value: 0.4 });
    expect(L.rate(complexes, "baseline", "valid").k).toBe(9);
    expect(L.rate(complexes, "box_padding_2", "top1").k).toBe(14);
    expect(L.rate(complexes, "box_padding_2", "valid").k).toBe(13);
  });

  it("counts gained and lost complexes", () => {
    expect(L.discordance(complexes, "box_padding_2", "top1")).toEqual({ gain: 5, loss: 1 });
    expect(L.discordance(complexes, "ligand_dimorphite", "top1")).toEqual({ gain: 0, loss: 0 });
  });

  it("orders the chart with the baseline first, then by effect", () => {
    const order = L.chartOrder(complexes, arms);
    expect(order.map((r) => r.arm.key)).toEqual(["baseline", "box_padding_2", "ligand_dimorphite"]);
    expect(order[0].delta).toBeNull();
    expect(order[1].delta).toBeCloseTo(0.16, 10);
    expect(order[2].delta).toBe(0);
  });

  it("cross-tabulates outcome classes", () => {
    const t = L.transitionCounts(complexes, "box_padding_2");
    expect(t.scoring.valid).toBe(5);
    expect(t.valid.mode).toBe(1);
    expect(t.valid.valid).toBe(8);
    expect(t.invalid.invalid).toBe(1);
    const total = L.CLASSES.reduce((s, a) => s + L.CLASSES.reduce((r, b) => r + t[a][b], 0), 0);
    expect(total).toBe(25);
  });
});

describe("filterComplexes", () => {
  const { complexes } = makeData();
  const base = { arm: "box_padding_2", metric: "top1", sort: { key: "id", dir: 1 } };
  const ids = (rows) => rows.map((c) => c.id);

  it("searches case-insensitively", () => {
    expect(ids(L.filterComplexes(complexes, { ...base, search: " c1" }))).toHaveLength(10);
    expect(ids(L.filterComplexes(complexes, { ...base, search: "c07" }))).toEqual(["C07_LIG"]);
  });

  it("filters by cofactor, outcome classes and changed success", () => {
    expect(L.filterComplexes(complexes, { ...base, cofactorOnly: true })).toHaveLength(5);
    expect(L.filterComplexes(complexes, { ...base, klass: "valid" })).toHaveLength(13);
    expect(L.filterComplexes(complexes, { ...base, baseClass: "scoring", klass: "valid" })).toHaveLength(5);
    expect(ids(L.filterComplexes(complexes, { ...base, changedOnly: true }))).toEqual([
      "C00_LIG", "C10_LIG", "C11_LIG", "C12_LIG", "C13_LIG", "C14_LIG",
    ]);
  });

  it("sorts by any column without mutating the input", () => {
    const before = ids(complexes);
    const rows = L.filterComplexes(complexes, { ...base, sort: { key: "heavy", dir: -1 } });
    expect(rows[0].id).toBe("C24_LIG");
    expect(rows.at(-1).id).toBe("C00_LIG");
    expect(ids(complexes)).toEqual(before);
  });
});

describe("paginate", () => {
  it.each([
    [25, 1, { page: 1, pages: 3, start: 0, end: 10 }],
    [25, 3, { page: 3, pages: 3, start: 20, end: 25 }],
    [25, Infinity, { page: 3, pages: 3, start: 20, end: 25 }],
    [25, 0, { page: 1, pages: 3, start: 0, end: 10 }],
    [0, 4, { page: 1, pages: 1, start: 0, end: 0 }],
    [10, 2, { page: 1, pages: 1, start: 0, end: 10 }],
  ])("%i rows, page %s", (total, page, expected) => {
    expect(L.paginate(total, page)).toEqual(expected);
  });

  it("finds the page holding a row", () => {
    expect([0, 9, 10, 24].map((i) => L.pageOf(i))).toEqual([1, 1, 2, 3]);
  });
});

describe("formatting", () => {
  it("formats percentages, signed changes and p-values", () => {
    expect(L.pct(0.512)).toBe("51.2%");
    expect(L.signed(13.1)).toBe("+13.1");
    expect(L.signed(-1.4)).toBe("−1.4"); // typographic minus
    expect(L.fmtP(0.0004)).toBe("< 0.001");
    expect(L.fmtP(0.0312)).toBe("0.031");
  });

  it("formats distances, with a dash when a complex has none", () => {
    expect(L.fmtA(1.234)).toBe("1.23");
    expect(L.fmtA(undefined)).toBe("–");
    expect(L.fmtA(null)).toBe("–");
  });
});

describe("pairedDiff", () => {
  it("matches the Python paired_difference_interval", () => {
    const [lo, hi] = L.pairedDiff(180, 18, 48, 428);
    expect(100 * lo).toBeCloseTo(3.34, 1);
    expect(100 * hi).toBeCloseTo(10.63, 1);
  });

  it("is symmetric around zero when nothing changes, and rejects bad counts", () => {
    const [lo, hi] = L.pairedDiff(100, 0, 0, 428);
    expect(lo).toBeCloseTo(-hi, 10);
    expect(() => L.pairedDiff(400, 20, 20, 428)).toThrow(RangeError);
  });
});

describe("agreement filter", () => {
  const base = { rmsd: 1, best: 1, valid: true, shift: 0 };
  const complexes = [
    { id: "A", cofactors: [], agree: 5, arms: { baseline: base } },
    { id: "B", cofactors: [], agree: 0, arms: { baseline: base } },
    { id: "C", cofactors: [], agree: null, arms: { baseline: base } },
  ];
  const f = (agree) => L.filterComplexes(complexes, { agree, arm: "baseline", metric: "top1" }).map((c) => c.id);

  it("keeps every complex when no level is chosen", () => {
    expect(f("")).toEqual(["A", "B", "C"]);
    expect(f(undefined)).toEqual(["A", "B", "C"]);
  });

  it("keeps only complexes at the chosen level, including 0", () => {
    expect(f("5")).toEqual(["A"]);
    expect(f("0")).toEqual(["B"]);
  });
});
