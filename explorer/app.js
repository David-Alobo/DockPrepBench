"use strict";

const DATA = window.EXPLORER_DATA;
const ARMS = DATA.arms;
const BASE = "baseline";
const CLASSES = ["valid", "invalid", "scoring", "site", "mode"];
const CLASS_LABEL = {
  valid: "Valid success",
  invalid: "Correct, not PB-valid",
  scoring: "Scoring miss",
  site: "Wrong site",
  mode: "Wrong binding mode",
};
const CLASS_SHORT = { valid: "Valid", invalid: "Not valid", scoring: "Scoring", site: "Site", mode: "Mode" };
const SITE_SHIFT = 4; // Å, heuristic split between wrong site and wrong binding mode
const PAGE_SIZE = 10;

const state = {
  tab: "overview",
  arm: "box_padding_2",
  metric: "top1",
  search: "",
  klass: "",
  baseClass: "",
  changedOnly: false,
  cofactorOnly: false,
  sort: { key: "id", dir: 1 },
  page: 1,
  selected: "5S8I_2LY",
  pose: 1,
  showCrystal: true,
  showPocket: true,
};

const $ = (id) => document.getElementById(id);
const css = (name) => getComputedStyle(document.documentElement).getPropertyValue(name).trim();
const armLabel = (key) => ARMS.find((a) => a.key === key).label;
const byId = new Map(DATA.complexes.map((c) => [c.id, c]));

// ---------------------------------------------------------------- statistics

function classify(o) {
  if (o.rmsd < 2) return o.valid ? "valid" : "invalid";
  if (o.best < 2) return "scoring";
  return o.shift > SITE_SHIFT ? "site" : "mode";
}

function succeeds(o, metric) {
  return o.rmsd < 2 && (metric === "top1" || o.valid);
}

function wilson(k, n, z = 1.959963984540054) {
  const p = k / n;
  const centre = (p + (z * z) / (2 * n)) / (1 + (z * z) / n);
  const half = (z * Math.sqrt((p * (1 - p)) / n + (z * z) / (4 * n * n))) / (1 + (z * z) / n);
  return [Math.max(0, centre - half), Math.min(1, centre + half)];
}

function logChoose(n, k) {
  let s = 0;
  for (let i = 1; i <= k; i++) s += Math.log(n - k + i) - Math.log(i);
  return s;
}

/** Exact two-sided McNemar p-value: binomial test on the discordant pairs. */
function mcnemar(gain, loss) {
  const n = gain + loss;
  if (n === 0) return 1;
  const k = Math.min(gain, loss);
  let tail = 0;
  for (let i = 0; i <= k; i++) tail += Math.exp(logChoose(n, i) - n * Math.LN2);
  return Math.min(1, 2 * tail);
}

function rate(arm, metric) {
  const n = DATA.complexes.length;
  const k = DATA.complexes.filter((c) => succeeds(c.arms[arm], metric)).length;
  const [lo, hi] = wilson(k, n);
  return { k, n, value: k / n, lo, hi };
}

function discordance() {
  let gain = 0, loss = 0;
  for (const c of DATA.complexes) {
    const a = succeeds(c.arms[BASE], state.metric);
    const b = succeeds(c.arms[state.arm], state.metric);
    if (b && !a) gain++;
    if (a && !b) loss++;
  }
  return { gain, loss };
}

const pct = (x) => `${(100 * x).toFixed(1)}%`;
const signed = (x, digits = 1) => `${x >= 0 ? "+" : "−"}${Math.abs(x).toFixed(digits)}`;
const fmtP = (p) => (p < 0.001 ? "< 0.001" : p.toPrecision(2));

// ---------------------------------------------------------------- glyphs

const glyph = (k, extra = "") => `<span class="g ${k}${extra}" aria-hidden="true"></span>`;
const chip = (k) => `<span class="chip">${glyph(k)}${CLASS_LABEL[k]}</span>`;

function strip(c) {
  const cells = ARMS.map((arm) => {
    const o = c.arms[arm.key];
    const k = classify(o);
    return `<span class="g ${k}${arm.key === state.arm ? " sel" : ""}" title="${arm.label}: ${CLASS_LABEL[k]}, top-1 ${o.rmsd.toFixed(2)} Å"></span>`;
  });
  const summary = ARMS.map((a) => `${a.label} ${CLASS_LABEL[classify(c.arms[a.key])]}`).join("; ");
  return `<span class="strip" role="img" aria-label="${summary}">${cells.join("")}</span>`;
}

// ---------------------------------------------------------------- overview: tiles

function drawTiles() {
  const base = rate(BASE, state.metric), arm = rate(state.arm, state.metric);
  const { gain, loss } = discordance();
  $("t-base").textContent = pct(base.value);
  $("t-base-ci").textContent = `95% CI ${pct(base.lo)}–${pct(base.hi)} · ${base.k}/${base.n}`;
  $("t-arm-k").textContent = `${armLabel(state.arm)} success`;
  $("t-arm").textContent = `${pct(arm.value)}`;
  $("t-arm-ci").textContent = `${signed(100 * (arm.value - base.value))} pts vs baseline · CI ${pct(arm.lo)}–${pct(arm.hi)}`;
  $("t-changed").textContent = `${gain + loss}`;
  $("t-changed-s").innerHTML = `<span class="delta up">↑ ${gain} gained</span> · <span class="delta down">↓ ${loss} lost</span>`;
  $("t-p").textContent = `p ${fmtP(mcnemar(gain, loss)).replace(/^(?!<)/, "= ")}`;
  $("t-p-s").textContent = `uncorrected; ${ARMS.length - 1} arms are compared with the baseline`;
}

// ---------------------------------------------------------------- overview: chart

function chartOrder() {
  const base = rate(BASE, "top1").value;
  const others = ARMS.filter((a) => a.key !== BASE)
    .map((a) => ({ arm: a, delta: rate(a.key, "top1").value - base }))
    .sort((a, b) => b.delta - a.delta);
  return [{ arm: ARMS.find((a) => a.key === BASE), delta: null }, ...others];
}

function drawChart() {
  const svg = $("chart");
  const W = 620, left = 126, right = 76, top = 22, rowH = 36, gap = 10;
  const order = chartOrder();
  const H = top + order.length * rowH + gap + 6;
  const x = (v) => left + v * (W - left - right);
  const rowY = (i) => top + i * rowH + (i > 0 ? gap : 0);
  svg.setAttribute("viewBox", `0 0 ${W} ${H}`);
  const parts = [];
  for (let t = 0; t <= 100; t += 20) {
    parts.push(`<line class="grid" x1="${x(t / 100)}" x2="${x(t / 100)}" y1="${top - 4}" y2="${H - 4}" />`);
    parts.push(`<text x="${x(t / 100)}" y="${top - 9}" text-anchor="middle">${t}%</text>`);
  }
  parts.push(`<text x="${W}" y="${top - 9}" text-anchor="end">Δ top-1</text>`);
  parts.push(`<line class="sep" x1="0" x2="${W}" y1="${rowY(1) - gap / 2 - 1}" y2="${rowY(1) - gap / 2 - 1}" />`);
  const ref = rate(BASE, "top1").value;
  parts.push(`<line class="ref" x1="${x(ref)}" x2="${x(ref)}" y1="${top - 4}" y2="${H - 4}" />`);
  const metrics = [
    ["top1", css("--metric-1"), -5, "Top-1 RMSD < 2 Å"],
    ["valid", css("--metric-2"), 5, "PB-valid success"],
  ];
  order.forEach(({ arm, delta }, i) => {
    const y0 = rowY(i), cy = y0 + rowH / 2;
    const isBase = arm.key === BASE, isSel = arm.key === state.arm;
    if (!isBase) {
      parts.push(`<rect class="row-hit" data-arm="${arm.key}" x="0" y="${y0}" width="${W}" height="${rowH}" rx="6" fill="${isSel ? css("--hover") : "transparent"}" style="cursor:pointer" />`);
    }
    const name = isBase ? "Baseline (ref.)" : arm.label;
    parts.push(`<text class="arm-label" x="6" y="${cy + 4}" font-weight="${isBase || isSel ? 600 : 400}" pointer-events="none">${name}</text>`);
    for (const [metric, colour, dy, label] of metrics) {
      const r = rate(arm.key, metric);
      const tip = `${name}<br>${label}: ${pct(r.value)} [${pct(r.lo)}, ${pct(r.hi)}]<br>${r.k} of ${r.n} complexes`;
      parts.push(
        `<g data-tip="${tip}" data-arm="${arm.key}">` +
          `<rect x="${x(r.lo) - 4}" y="${cy + dy - 6}" width="${x(r.hi) - x(r.lo) + 8}" height="12" fill="transparent" />` +
          `<line x1="${x(r.lo)}" x2="${x(r.hi)}" y1="${cy + dy}" y2="${cy + dy}" stroke="${colour}" stroke-width="2" stroke-linecap="round" pointer-events="none" />` +
          `<circle cx="${x(r.value)}" cy="${cy + dy}" r="4.5" fill="${colour}" stroke="${css("--surface")}" stroke-width="2" pointer-events="none" /></g>`,
      );
    }
    const dText = delta === null ? "ref" : `${signed(100 * delta)} pts`;
    parts.push(`<text class="delta-label" x="${W - 4}" y="${cy + 4}" text-anchor="end" pointer-events="none">${dText}</text>`);
  });
  svg.innerHTML = parts.join("");
  svg.querySelectorAll("[data-tip]").forEach((g) => {
    g.addEventListener("mousemove", (e) => showTip(e, g.dataset.tip));
    g.addEventListener("mouseleave", hideTip);
  });
  svg.querySelectorAll("[data-arm]").forEach((el) =>
    el.addEventListener("click", () => {
      if (el.dataset.arm !== BASE) setArm(el.dataset.arm);
    }),
  );
}

function showTip(e, html) {
  const tip = $("tip");
  tip.innerHTML = html;
  tip.hidden = false;
  tip.style.left = `${Math.min(e.clientX + 14, window.innerWidth - tip.offsetWidth - 8)}px`;
  tip.style.top = `${e.clientY + 14}px`;
}

function hideTip() {
  $("tip").hidden = true;
}

// ---------------------------------------------------------------- overview: transition matrix

function drawMatrix() {
  const counts = Object.fromEntries(CLASSES.map((a) => [a, Object.fromEntries(CLASSES.map((b) => [b, 0]))]));
  for (const c of DATA.complexes) counts[classify(c.arms[BASE])][classify(c.arms[state.arm])]++;
  const max = Math.max(...CLASSES.flatMap((a) => CLASSES.map((b) => counts[a][b])));
  const head = `<tr><th class="row">Baseline ↓ · ${armLabel(state.arm)} →</th>${CLASSES.map(
    (b) => `<th title="${CLASS_LABEL[b]}">${glyph(b)} ${CLASS_SHORT[b]}</th>`,
  ).join("")}</tr>`;
  const body = CLASSES.map((a) => {
    const cells = CLASSES.map((b) => {
      const n = counts[a][b];
      if (!n) return `<td class="zero${a === b ? " diag" : ""}">0</td>`;
      const share = Math.round(15 + 70 * Math.sqrt(n / max));
      const ink = share > 55 ? "#ffffff" : "var(--fg)";
      const label = `${n} complexes: ${CLASS_LABEL[a]} at baseline, ${CLASS_LABEL[b]} with ${armLabel(state.arm)}`;
      return (
        `<td data-cell="${a}|${b}" tabindex="0" title="${label}" aria-label="${label}" class="${a === b ? "diag" : ""}" ` +
        `style="background: color-mix(in srgb, var(--heat-1) ${share}%, var(--heat-0)); color: ${ink}">${n}</td>`
      );
    }).join("");
    return `<tr><th class="row" title="${CLASS_LABEL[a]}">${glyph(a)} ${CLASS_LABEL[a]}</th>${cells}</tr>`;
  }).join("");
  $("matrix").innerHTML = `<thead>${head}</thead><tbody>${body}</tbody>`;
  $("matrix-axes").textContent = `${DATA.complexes.length} complexes`;
}

function openCell(cell) {
  const [a, b] = cell.split("|");
  state.baseClass = a;
  state.klass = b;
  state.changedOnly = false;
  state.search = "";
  $("class-filter").value = b;
  $("changed-only").checked = false;
  $("search").value = "";
  state.page = 1;
  const first = filteredRows()[0];
  if (first) state.selected = first.id;
  showTab("complexes");
}

// ---------------------------------------------------------------- complexes: table

function filteredRows() {
  const q = state.search.trim().toUpperCase();
  const rows = DATA.complexes.filter((c) => {
    if (q && !c.id.includes(q)) return false;
    if (state.cofactorOnly && !c.cofactors.length) return false;
    if (state.klass && classify(c.arms[state.arm]) !== state.klass) return false;
    if (state.baseClass && classify(c.arms[BASE]) !== state.baseClass) return false;
    if (state.changedOnly) {
      return succeeds(c.arms[BASE], state.metric) !== succeeds(c.arms[state.arm], state.metric);
    }
    return true;
  });
  const { key, dir } = state.sort;
  return rows.sort((a, b) => (a[key] < b[key] ? -dir : a[key] > b[key] ? dir : 0));
}

function drawTable() {
  const rows = filteredRows();
  const pages = Math.max(1, Math.ceil(rows.length / PAGE_SIZE));
  state.page = Math.min(Math.max(1, state.page), pages);
  const start = (state.page - 1) * PAGE_SIZE;
  const shown = rows.slice(start, start + PAGE_SIZE);
  const filler = Array.from({ length: PAGE_SIZE - shown.length }, () => `<tr aria-hidden="true"><td colspan="6"></td></tr>`);
  $("rows").innerHTML = shown.length
    ? shown
        .map((c) => {
          const a = succeeds(c.arms[BASE], state.metric);
          const b = succeeds(c.arms[state.arm], state.metric);
          const change = a === b ? `<span class="delta same">same</span>`
            : b ? `<span class="delta up">↑ gained</span>` : `<span class="delta down">↓ lost</span>`;
          return (
            `<tr tabindex="0" data-id="${c.id}" aria-selected="${c.id === state.selected}">` +
            `<td class="mono">${c.id}${DATA.real[c.id] ? " ◆" : ""}</td><td class="num">${c.heavy}</td><td class="num">${c.rotb}</td>` +
            `<td>${strip(c)}</td><td>${chip(classify(c.arms[state.arm]))}</td><td>${change}</td></tr>`
          );
        })
        .concat(filler)
        .join("")
    : `<tr><td colspan="6" class="empty-row">No complexes match these filters.</td></tr>` + filler.slice(1).join("");

  const end = Math.min(start + PAGE_SIZE, rows.length);
  $("count").textContent = rows.length
    ? `${start + 1}–${end} of ${rows.length}${rows.length < DATA.complexes.length ? ` (filtered from ${DATA.complexes.length})` : ""}`
    : "none match";
  $("page-label").textContent = `${state.page} / ${pages}`;
  $("page-first").disabled = $("page-prev").disabled = state.page === 1;
  $("page-next").disabled = $("page-last").disabled = state.page === pages;
  document.querySelectorAll(".list th button[data-sort]").forEach((b) => {
    b.dataset.active = b.dataset.sort === state.sort.key ? String(state.sort.dir) : "";
  });
  $("arm-col").textContent = armLabel(state.arm);
  $("strip-head").title = `Baseline, then arms: ${ARMS.slice(1).map((a) => a.label).join(", ")}`;

  const af = $("active-filter");
  af.hidden = !state.baseClass;
  if (state.baseClass) {
    af.innerHTML =
      `<span class="muted">Showing</span> ${glyph(state.baseClass)} ${CLASS_LABEL[state.baseClass]} at baseline` +
      (state.klass ? ` → ${glyph(state.klass)} ${CLASS_LABEL[state.klass]} with ${armLabel(state.arm)}` : "") +
      ` <button type="button" id="clear-cell">Clear</button>`;
    $("clear-cell").addEventListener("click", () => {
      state.baseClass = "";
      state.klass = "";
      $("class-filter").value = "";
      refilter();
    });
  }
}

/** Re-render from page 1 after anything that changes which rows match. */
function refilter() {
  state.page = 1;
  drawTable();
}

function revealSelected() {
  const i = filteredRows().findIndex((c) => c.id === state.selected);
  if (i >= 0) state.page = Math.floor(i / PAGE_SIZE) + 1;
}

// ---------------------------------------------------------------- complexes: detail and 3D viewer

let viewer = null;

function drawDetail() {
  const c = byId.get(state.selected);
  if (!c) return;
  $("d-id").textContent = c.id;
  $("d-depict").innerHTML = c.svg;
  $("d-props").innerHTML =
    `<span>Heavy atoms <b>${c.heavy}</b></span><span>Rot. bonds <b>${c.rotb}</b></span>` +
    `<span>MW <b>${c.mw}</b></span><span>Cofactors <b>${c.cofactors.length ? c.cofactors.join(", ") : "none"}</b></span>`;
  $("d-arms").innerHTML = ARMS.map((arm) => {
    const o = c.arms[arm.key];
    return (
      `<div class="arm-cell${arm.key === state.arm ? " sel" : ""}"><div class="name">${arm.label}${o.real ? " ◆" : ""}</div>` +
      `<div class="val">${glyph(classify(o))}<span>${o.rmsd.toFixed(2)} Å</span>` +
      `<span class="best" title="Best RMSD among all returned poses">best ${o.best.toFixed(2)}</span></div></div>`
    );
  }).join("");

  const real = DATA.real[c.id];
  $("viewer-bar").hidden = !real;
  if (!real) {
    viewer = null;
    $("viewer").innerHTML =
      `<div class="viewer-empty">Docked poses for this complex will appear here after the full benchmark run.<br>Complexes marked ◆ have real poses now.</div>`;
    $("d-poses").innerHTML = `<span class="count">No poses saved yet.</span>`;
    return;
  }
  $("d-poses").innerHTML = real.poses
    .map((p) => {
      const ok = p.rmsd < 2 ? ` <span class="ok">✓</span>` : "";
      return (
        `<button type="button" data-rank="${p.rank}" aria-pressed="${p.rank === state.pose}" aria-label="Pose ${p.rank}, score ${p.score}, RMSD ${p.rmsd} Å">` +
        `<span class="rank">#${p.rank}</span><span>${p.score.toFixed(1)}</span><span>${p.rmsd.toFixed(2)}${ok}</span></button>`
      );
    })
    .join("");
  if (state.tab === "complexes") drawViewer(real);
}

function drawViewer(real) {
  const el = $("viewer");
  if (typeof $3Dmol === "undefined") {
    el.innerHTML = `<div class="viewer-empty">The 3D viewer library did not load. Check your connection and reload.</div>`;
    return;
  }
  if (!viewer) {
    el.innerHTML = "";
    viewer = $3Dmol.createViewer(el, { backgroundColor: css("--viewer-bg"), antialias: true });
  }
  viewer.setBackgroundColor(css("--viewer-bg"));
  viewer.clear();
  if (state.showPocket) {
    const pocket = viewer.addModel(real.pocket, "pdb");
    viewer.setStyle({ model: pocket }, { stick: { radius: 0.08, colorscheme: { prop: "elem", map: elementColours("#9aa3b2") } } });
  }
  const crystal = viewer.addModel(real.crystal, "sdf");
  viewer.setStyle({ model: crystal }, state.showCrystal ? { stick: { radius: 0.2, color: css("--faint"), opacity: 0.85 } } : {});
  const pose = real.poses.find((p) => p.rank === state.pose) || real.poses[0];
  const docked = viewer.addModel(pose.block, "sdf");
  viewer.setStyle({ model: docked }, { stick: { radius: 0.24, colorscheme: { prop: "elem", map: elementColours(css("--metric-1")) } } });
  viewer.resize();
  viewer.zoomTo({ model: crystal });
  viewer.zoom(0.8);
  viewer.render();
}

function elementColours(carbon) {
  return { C: carbon, N: "#3b5bdb", O: "#e03131", S: "#e0a800", F: "#2f9e44", Cl: "#2f9e44", Br: "#a5541f", P: "#e8590c", H: "#dee2e6" };
}

function select(id) {
  if (!byId.has(id)) return;
  state.selected = id;
  state.pose = 1;
  try {
    history.replaceState(null, "", `#${id}`);
  } catch {
    /* sandboxed frames may refuse; the copy-link button still works */
  }
  drawTable();
  drawDetail();
}

async function copyLink() {
  const base = DATA.shareUrl || `${location.origin}${location.pathname}`;
  const link = `${base}#${state.selected}`;
  const button = $("copy-link");
  try {
    await navigator.clipboard.writeText(link);
    button.textContent = "Copied";
  } catch {
    button.textContent = link;
  }
  setTimeout(() => { button.textContent = "Copy link"; }, 2500);
}

// ---------------------------------------------------------------- wiring

function showTab(tab) {
  state.tab = tab;
  for (const name of ["overview", "complexes"]) {
    $(`tab-${name}`).setAttribute("aria-selected", String(name === tab));
    $(`panel-${name}`).hidden = name !== tab;
  }
  if (tab === "complexes") {
    revealSelected();
    drawTable();
    drawDetail(); // the 3D viewer needs a visible container to size itself
  }
}

function setArm(key) {
  state.arm = key;
  state.page = 1;
  $("arm-select").value = key;
  renderAll();
}

function renderAll() {
  drawTiles();
  drawChart();
  drawMatrix();
  drawTable();
  drawDetail();
}

function init() {
  $("meta").textContent = `${DATA.dataset} · ${DATA.complexes.length} complexes · ${DATA.engine}`;
  $("badge").hidden = !DATA.demo;
  $("key").innerHTML = CLASSES.map((k) => `<span>${glyph(k)} ${CLASS_LABEL[k]}</span>`).join("");
  $("footer").textContent =
    "Outcome classes follow the top-ranked pose: valid success (RMSD < 2 Å and PB-valid); correct but not PB-valid; " +
    `scoring miss (a pose under 2 Å exists below rank 1); wrong site (centroid more than ${SITE_SHIFT} Å away, a heuristic); ` +
    "wrong binding mode (right site, wrong pose).";
  $("arm-select").innerHTML = ARMS.filter((a) => a.key !== BASE)
    .map((a) => `<option value="${a.key}">${a.label}: ${a.detail}</option>`)
    .join("");
  $("arm-select").value = state.arm;

  $("tab-overview").addEventListener("click", () => showTab("overview"));
  $("tab-complexes").addEventListener("click", () => showTab("complexes"));
  $("arm-select").addEventListener("change", (e) => setArm(e.target.value));
  $("metric-select").addEventListener("change", (e) => { state.metric = e.target.value; state.page = 1; renderAll(); });
  $("search").addEventListener("input", (e) => { state.search = e.target.value; refilter(); });
  $("class-filter").addEventListener("change", (e) => { state.klass = e.target.value; refilter(); });
  $("changed-only").addEventListener("change", (e) => { state.changedOnly = e.target.checked; refilter(); });
  $("cofactor-only").addEventListener("change", (e) => { state.cofactorOnly = e.target.checked; refilter(); });
  document.querySelectorAll(".list th button[data-sort]").forEach((b) =>
    b.addEventListener("click", () => {
      const key = b.dataset.sort;
      state.sort = { key, dir: state.sort.key === key ? -state.sort.dir : 1 };
      refilter();
    }),
  );
  const goTo = (page) => { state.page = page; drawTable(); };
  $("page-first").addEventListener("click", () => goTo(1));
  $("page-prev").addEventListener("click", () => goTo(state.page - 1));
  $("page-next").addEventListener("click", () => goTo(state.page + 1));
  $("page-last").addEventListener("click", () => goTo(Infinity)); // clamped to the last page

  $("matrix").addEventListener("click", (e) => {
    const td = e.target.closest("td[data-cell]");
    if (td) openCell(td.dataset.cell);
  });
  $("matrix").addEventListener("keydown", (e) => {
    const td = e.target.closest("td[data-cell]");
    if (td && (e.key === "Enter" || e.key === " ")) { e.preventDefault(); openCell(td.dataset.cell); }
  });
  $("rows").addEventListener("click", (e) => {
    const tr = e.target.closest("tr[data-id]");
    if (tr) select(tr.dataset.id);
  });
  $("rows").addEventListener("keydown", (e) => {
    const tr = e.target.closest("tr[data-id]");
    if (tr && (e.key === "Enter" || e.key === " ")) { e.preventDefault(); select(tr.dataset.id); }
  });
  $("d-poses").addEventListener("click", (e) => {
    const b = e.target.closest("button[data-rank]");
    if (!b) return;
    state.pose = Number(b.dataset.rank);
    drawDetail();
  });
  $("t-crystal").addEventListener("change", (e) => { state.showCrystal = e.target.checked; drawDetail(); });
  $("t-pocket").addEventListener("change", (e) => { state.showPocket = e.target.checked; drawDetail(); });
  $("copy-link").addEventListener("click", copyLink);

  const redraw = () => { drawChart(); drawMatrix(); drawDetail(); };
  window.matchMedia("(prefers-color-scheme: dark)").addEventListener("change", redraw);
  new MutationObserver(redraw).observe(document.documentElement, { attributes: true, attributeFilter: ["data-theme"] });
  window.addEventListener("resize", () => viewer && viewer.resize());

  renderAll();
  // A link to a complex (e.g. #7KLX_WOV) opens it on the Complexes tab.
  const fromLink = decodeURIComponent(location.hash.slice(1)).toUpperCase();
  if (byId.has(fromLink)) {
    state.selected = fromLink;
    showTab("complexes");
  }
}

init();
