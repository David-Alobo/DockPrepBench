/**
 * Pure logic for the results explorer: outcome classes, statistics, filtering and paging.
 *
 * No DOM access, so it can be unit-tested. Loaded as a plain <script> in the page (it sets
 * `window.ExplorerLib`) and as a CommonJS module in tests.
 */
(function (root, factory) {
  const lib = factory();
  if (typeof module === "object" && module.exports) module.exports = lib;
  else root.ExplorerLib = lib;
})(typeof globalThis !== "undefined" ? globalThis : this, function () {
  "use strict";

  const BASE = "baseline";
  const CLASSES = ["valid", "invalid", "scoring", "site", "mode", "failed"];
  const CLASS_LABEL = {
    valid: "Valid success",
    invalid: "Correct, not PB-valid",
    scoring: "Scoring miss",
    site: "Wrong site",
    mode: "Wrong binding mode",
    failed: "Failed to prepare or dock",
  };
  const CLASS_SHORT = { valid: "Valid", invalid: "Not valid", scoring: "Scoring", site: "Site", mode: "Mode", failed: "Failed" };
  const SUCCESS_RMSD = 2; // Å
  const SITE_SHIFT = 4; // Å, heuristic split between wrong site and wrong binding mode
  const PAGE_SIZE = 10;
  const Z_95 = 1.959963984540054;

  /**
   * Classify a top-1 outcome: failed (preparation or docking raised an error), valid,
   * invalid (correct but not PB-valid), scoring miss, wrong site or wrong binding mode.
   * The first matching rule wins.
   * @param {{rmsd?: number, best?: number, valid?: boolean, shift?: number, error?: string}} o
   * @param {number} [siteShift] Centroid shift (Å) above which a miss is a wrong site.
   * @returns {string} One of CLASSES.
   */
  function classify(o, siteShift = SITE_SHIFT) {
    if (o.error) return "failed";
    if (o.rmsd < SUCCESS_RMSD) return o.valid ? "valid" : "invalid";
    if (o.best < SUCCESS_RMSD) return "scoring";
    return o.shift > siteShift ? "site" : "mode";
  }

  /**
   * Whether an outcome counts as a success under a metric.
   * @param {{rmsd?: number, valid?: boolean, error?: string}} o A failed complex never succeeds.
   * @param {"top1"|"valid"} metric "top1": RMSD < 2 Å; "valid": also PB-valid.
   * @returns {boolean}
   */
  function succeeds(o, metric) {
    return !o.error && o.rmsd < SUCCESS_RMSD && (metric === "top1" || o.valid);
  }

  /**
   * Wilson score interval for a binomial proportion (matches analysis.stats.wilson_interval).
   * @param {number} k Successes.
   * @param {number} n Trials (> 0).
   * @param {number} [z] Normal quantile; default gives a 95% interval.
   * @returns {[number, number]} Lower and upper bounds.
   */
  function wilson(k, n, z = Z_95) {
    if (!(n > 0) || k < 0 || k > n) throw new RangeError(`Need 0 <= k <= n and n > 0, got k=${k}, n=${n}`);
    const p = k / n;
    const centre = (p + (z * z) / (2 * n)) / (1 + (z * z) / n);
    const half = (z * Math.sqrt((p * (1 - p)) / n + (z * z) / (4 * n * n))) / (1 + (z * z) / n);
    return [Math.max(0, centre - half), Math.min(1, centre + half)];
  }

  /** Natural log of the binomial coefficient C(n, k). */
  function logChoose(n, k) {
    let s = 0;
    for (let i = 1; i <= k; i++) s += Math.log(n - k + i) - Math.log(i);
    return s;
  }

  /**
   * Exact two-sided McNemar p-value: a binomial test (p = 0.5) on the discordant pairs.
   * @param {number} gain Complexes only the arm solves.
   * @param {number} loss Complexes only the baseline solves.
   * @returns {number} p-value in [0, 1].
   */
  function mcnemar(gain, loss) {
    const n = gain + loss;
    if (n === 0) return 1;
    const k = Math.min(gain, loss);
    let tail = 0;
    for (let i = 0; i <= k; i++) tail += Math.exp(logChoose(n, i) - n * Math.LN2);
    return Math.min(1, 2 * tail);
  }

  /**
   * Success rate of one arm with its Wilson interval.
   * @param {Array} complexes Complexes with `arms[arm]` outcomes.
   * @param {string} arm Arm key.
   * @param {"top1"|"valid"} metric
   * @returns {{k: number, n: number, value: number, lo: number, hi: number}}
   */
  function rate(complexes, arm, metric) {
    const n = complexes.length;
    const k = complexes.filter((c) => succeeds(c.arms[arm], metric)).length;
    const [lo, hi] = wilson(k, n);
    return { k, n, value: k / n, lo, hi };
  }

  /**
   * Count complexes whose success differs between the baseline and an arm.
   * @returns {{gain: number, loss: number}} gain: only the arm succeeds; loss: only the baseline.
   */
  function discordance(complexes, arm, metric, base = BASE) {
    let gain = 0;
    let loss = 0;
    for (const c of complexes) {
      const a = succeeds(c.arms[base], metric);
      const b = succeeds(c.arms[arm], metric);
      if (b && !a) gain++;
      if (a && !b) loss++;
    }
    return { gain, loss };
  }

  /**
   * Chart row order: the baseline first, then the arms by descending change in top-1 success.
   * @returns {Array<{arm: object, delta: number|null}>} delta is null for the baseline.
   */
  function chartOrder(complexes, arms, base = BASE) {
    const ref = rate(complexes, base, "top1").value;
    const others = arms
      .filter((a) => a.key !== base)
      .map((a) => ({ arm: a, delta: rate(complexes, a.key, "top1").value - ref }))
      .sort((a, b) => b.delta - a.delta);
    return [{ arm: arms.find((a) => a.key === base), delta: null }, ...others];
  }

  /**
   * Cross-tabulate outcome classes from the baseline (outer key) to an arm (inner key).
   * @returns {Object<string, Object<string, number>>}
   */
  function transitionCounts(complexes, arm, base = BASE) {
    const counts = Object.fromEntries(CLASSES.map((a) => [a, Object.fromEntries(CLASSES.map((b) => [b, 0]))]));
    for (const c of complexes) counts[classify(c.arms[base])][classify(c.arms[arm])]++;
    return counts;
  }

  /**
   * Filter and sort complexes for the table.
   * @param {Array} complexes
   * @param {object} f Filters: search, cofactorOnly, klass (arm class), baseClass,
   *   changedOnly, arm, metric, sort {key, dir}.
   * @returns {Array} A new, sorted array.
   */
  function filterComplexes(complexes, f, base = BASE) {
    const q = (f.search || "").trim().toUpperCase();
    const rows = complexes.filter((c) => {
      if (q && !c.id.toUpperCase().includes(q)) return false;
      if (f.cofactorOnly && !c.cofactors.length) return false;
      if (f.klass && classify(c.arms[f.arm]) !== f.klass) return false;
      if (f.baseClass && classify(c.arms[base]) !== f.baseClass) return false;
      if (f.changedOnly) return succeeds(c.arms[base], f.metric) !== succeeds(c.arms[f.arm], f.metric);
      return true;
    });
    const { key, dir } = f.sort || { key: "id", dir: 1 };
    return rows.sort((a, b) => (a[key] < b[key] ? -dir : a[key] > b[key] ? dir : 0));
  }

  /**
   * Clamp a page number and compute the slice of rows it shows.
   * @param {number} total Number of rows.
   * @param {number} page Requested 1-based page (may be out of range, e.g. Infinity).
   * @param {number} [size]
   * @returns {{page: number, pages: number, start: number, end: number}} end is exclusive.
   */
  function paginate(total, page, size = PAGE_SIZE) {
    const pages = Math.max(1, Math.ceil(total / size));
    const clamped = Math.min(Math.max(1, page), pages);
    const start = (clamped - 1) * size;
    return { page: clamped, pages, start, end: Math.min(start + size, total) };
  }

  /** 1-based page that contains row `index`. */
  function pageOf(index, size = PAGE_SIZE) {
    return Math.floor(index / size) + 1;
  }

  const pct = (x) => `${(100 * x).toFixed(1)}%`;
  const signed = (x, digits = 1) => `${x >= 0 ? "+" : "−"}${Math.abs(x).toFixed(digits)}`;
  /** p-value as "< 0.001" or 2 significant figures. */
  const fmtP = (p) => (p < 0.001 ? "< 0.001" : p.toPrecision(2));
  /** A distance in Å to 2 decimals, or a dash when there is none (failed complexes). */
  const fmtA = (x) => (typeof x === "number" && Number.isFinite(x) ? x.toFixed(2) : "–");

  return {
    BASE,
    CLASSES,
    CLASS_LABEL,
    CLASS_SHORT,
    SITE_SHIFT,
    PAGE_SIZE,
    classify,
    succeeds,
    wilson,
    mcnemar,
    rate,
    discordance,
    chartOrder,
    transitionCounts,
    filterComplexes,
    paginate,
    pageOf,
    pct,
    signed,
    fmtP,
    fmtA,
  };
});
