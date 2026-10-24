"use strict";

// Page wiring for the results explorer. All logic that can be tested without a page lives
// in lib.js (window.ExplorerLib).
const L = window.ExplorerLib;
const DATA = window.EXPLORER_DATA;
const { BASE, CLASSES, CLASS_LABEL, CLASS_SHORT, SITE_SHIFT, PAGE_SIZE } = L;

// Built data holds several engines, each with its own arms, and per-complex outcomes in
// c.byEngine; demo and test data hold one engine, with outcomes already in c.arms.
// useEngine() points ARMS and every c.arms at one engine, so the rest of the page is unchanged.
const ENGINES = DATA.engines || null;
let ARMS = DATA.arms;
const defaultArm = (arms) =>
  arms.some((a) => a.key === "box_padding_2") ? "box_padding_2" : (arms.find((a) => a.key !== BASE) || arms[0]).key;

function useEngine(key) {
  const engine = ENGINES.find((e) => e.key === key) || ENGINES[0];
  ARMS = engine.arms;
  DATA.engine = engine.label;
  for (const c of DATA.complexes) c.arms = c.byEngine[engine.key];
  return engine.key;
}

const state = {
  tab: "overview",
  engine: ENGINES ? useEngine(new URLSearchParams(location.search).get("engine")) : null,
  arm: null, // set below, once ARMS is known

  metric: "top1",
  search: "",
  klass: "",
  baseClass: "",
  changedOnly: false,
  cofactorOnly: false,
  sort: { key: "id", dir: 1 },
  page: 1,
  selected: DATA.complexes.some((c) => c.id === "5S8I_2LY") ? "5S8I_2LY" : DATA.complexes[0].id,
  pose: 1,
  showCrystal: true,
  showPocket: true,
  showProtein: true,
  showSurface: true,
  closeUp: false,
};
state.arm = defaultArm(ARMS);

const $ = (id) => document.getElementById(id);
const css = (name) => getComputedStyle(document.documentElement).getPropertyValue(name).trim();
const armLabel = (key) => ARMS.find((a) => a.key === key).label;
const byId = new Map(DATA.complexes.map((c) => [c.id, c]));
const rate = (arm, metric) => L.rate(DATA.complexes, arm, metric);
const filters = () => ({ ...state });

// ---------------------------------------------------------------- glyphs

const glyph = (k, extra = "") => `<span class="g ${k}${extra}" aria-hidden="true"></span>`;
const chip = (k) => `<span class="chip" title="${CLASS_LABEL[k]}">${glyph(k)}${CLASS_LABEL[k]}</span>`;

function strip(c) {
  const cells = ARMS.map((arm) => {
    const o = c.arms[arm.key];
    const k = L.classify(o);
    const detail = o.error ? "" : `, top-1 ${L.fmtA(o.rmsd)} Å`;
    return `<span class="g ${k}${arm.key === state.arm ? " sel" : ""}" title="${arm.label}: ${CLASS_LABEL[k]}${detail}"></span>`;
  });
  const summary = ARMS.map((a) => `${a.label} ${CLASS_LABEL[L.classify(c.arms[a.key])]}`).join("; ");
  return `<span class="strip" role="img" aria-label="${summary}">${cells.join("")}</span>`;
}

// ---------------------------------------------------------------- overview: tiles

function drawTiles() {
  const base = rate(BASE, state.metric);
  const arm = rate(state.arm, state.metric);
  const { gain, loss } = L.discordance(DATA.complexes, state.arm, state.metric);
  const p = L.fmtP(L.mcnemar(gain, loss));
  $("t-base").textContent = L.pct(base.value);
  $("t-base-ci").textContent = `95% CI ${L.pct(base.lo)}–${L.pct(base.hi)} · ${base.k}/${base.n}`;
  $("t-arm-k").textContent = `${armLabel(state.arm)} success`;
  $("t-arm").textContent = L.pct(arm.value);
  $("t-arm-ci").textContent = `${L.signed(100 * (arm.value - base.value))} pts vs baseline · CI ${L.pct(arm.lo)}–${L.pct(arm.hi)}`;
  $("t-changed").textContent = `${gain + loss}`;
  $("t-changed-s").innerHTML = `<span class="delta up">↑ ${gain} gained</span> · <span class="delta down">↓ ${loss} lost</span>`;
  $("t-p").textContent = p.startsWith("<") ? `p ${p}` : `p = ${p}`;
  $("t-p-s").textContent = `uncorrected; ${ARMS.length - 1} arms are compared with the baseline`;
}

// ---------------------------------------------------------------- overview: chart

function drawChart() {
  const svg = $("chart");
  const W = 620, left = 126, right = 76, top = 22, rowH = 36, gap = 10;
  const order = L.chartOrder(DATA.complexes, ARMS);
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
      const tip = `${name}<br>${label}: ${L.pct(r.value)} [${L.pct(r.lo)}, ${L.pct(r.hi)}]<br>${r.k} of ${r.n} complexes`;
      parts.push(
        `<g data-tip="${tip}" data-arm="${arm.key}">` +
          `<rect x="${x(r.lo) - 4}" y="${cy + dy - 6}" width="${x(r.hi) - x(r.lo) + 8}" height="12" fill="transparent" />` +
          `<line x1="${x(r.lo)}" x2="${x(r.hi)}" y1="${cy + dy}" y2="${cy + dy}" stroke="${colour}" stroke-width="2" stroke-linecap="round" pointer-events="none" />` +
          `<circle cx="${x(r.value)}" cy="${cy + dy}" r="4.5" fill="${colour}" stroke="${css("--surface")}" stroke-width="2" pointer-events="none" /></g>`,
      );
    }
    const dText = delta === null ? "ref" : `${L.signed(100 * delta)} pts`;
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
  const counts = L.transitionCounts(DATA.complexes, state.arm);
  const max = Math.max(1, ...CLASSES.flatMap((a) => CLASSES.map((b) => counts[a][b])));
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
  Object.assign(state, { baseClass: a, klass: b, changedOnly: false, search: "", page: 1 });
  $("class-filter").value = b;
  $("changed-only").checked = false;
  $("search").value = "";
  const first = L.filterComplexes(DATA.complexes, filters())[0];
  if (first) state.selected = first.id;
  showTab("complexes");
}

// ---------------------------------------------------------------- complexes: table

function drawTable() {
  const rows = L.filterComplexes(DATA.complexes, filters());
  const { page, pages, start, end } = L.paginate(rows.length, state.page);
  state.page = page;
  const shown = rows.slice(start, end);
  const filler = Array.from({ length: PAGE_SIZE - shown.length }, () => `<tr aria-hidden="true"><td colspan="6"></td></tr>`);
  $("rows").innerHTML = shown.length
    ? shown
        .map((c) => {
          const a = L.succeeds(c.arms[BASE], state.metric);
          const b = L.succeeds(c.arms[state.arm], state.metric);
          const change = a === b ? `<span class="delta same">same</span>`
            : b ? `<span class="delta up">↑ gained</span>` : `<span class="delta down">↓ lost</span>`;
          return (
            `<tr tabindex="0" data-id="${c.id}" aria-selected="${c.id === state.selected}">` +
            `<td class="mono">${c.id}${DATA.real[c.id] ? " ◆" : ""}</td><td class="num">${c.heavy}</td><td class="num">${c.rotb}</td>` +
            `<td>${strip(c)}</td><td>${chip(L.classify(c.arms[state.arm]))}</td><td>${change}</td></tr>`
          );
        })
        .concat(filler)
        .join("")
    : `<tr><td colspan="6" class="empty-row">No complexes match these filters.</td></tr>` + filler.slice(1).join("");

  const total = DATA.complexes.length;
  $("count").textContent = rows.length
    ? `${start + 1}–${end} of ${rows.length}${rows.length < total ? ` (filtered from ${total})` : ""}`
    : "none match";
  $("page-label").textContent = `${page} / ${pages}`;
  $("page-first").disabled = $("page-prev").disabled = page === 1;
  $("page-next").disabled = $("page-last").disabled = page === pages;
  document.querySelectorAll(".list th button[data-sort]").forEach((b) => {
    b.dataset.active = b.dataset.sort === state.sort.key ? String(state.sort.dir) : "";
  });
  $("arm-col").textContent = armLabel(state.arm);
  $("strip-head").title = `Baseline, then arms: ${ARMS.filter((a) => a.key !== BASE).map((a) => a.label).join(", ")}`;

  const af = $("active-filter");
  af.hidden = !state.baseClass;
  if (state.baseClass) {
    af.innerHTML =
      `<span class="muted">Showing</span> ${glyph(state.baseClass)} ${CLASS_LABEL[state.baseClass]} at baseline` +
      (state.klass ? ` → ${glyph(state.klass)} ${CLASS_LABEL[state.klass]} with ${armLabel(state.arm)}` : "") +
      ` <button type="button" id="clear-cell">Clear</button>`;
    $("clear-cell").addEventListener("click", () => {
      Object.assign(state, { baseClass: "", klass: "" });
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
  const i = L.filterComplexes(DATA.complexes, filters()).findIndex((c) => c.id === state.selected);
  if (i >= 0) state.page = L.pageOf(i);
}

// ---------------------------------------------------------------- complexes: detail and 3D viewer

let viewer = null;

// Poses: demo data embeds them in DATA.real; built data keeps one script per complex under
// DATA.poseDir, which calls window.ExplorerPoses when it loads. null marks a missing file.
const poseCache = new Map(Object.entries(DATA.real || {}));
const posePending = new Set();

window.ExplorerPoses = (id, payload) => {
  posePending.delete(id);
  poseCache.set(id, payload);
  if (id === state.selected) drawDetail();
};

function requestPoses(id) {
  if (!DATA.poseDir || poseCache.has(id) || posePending.has(id)) return;
  posePending.add(id);
  const script = document.createElement("script");
  script.src = `${DATA.poseDir}/${encodeURIComponent(id)}.js`;
  script.onerror = () => window.ExplorerPoses(id, null);
  document.head.appendChild(script);
}

/** The selected arm's poses from a payload (demo payloads hold one list for the baseline). */
function armPoses(payload) {
  if (!payload) return [];
  if (Array.isArray(payload.poses)) return payload.poses;
  const byArm = ENGINES ? payload.poses[state.engine] || {} : payload.poses;
  return byArm[state.arm] || [];
}

function showNoPoses(message, list) {
  viewer = null;
  $("viewer-bar").hidden = true;
  $("viewer-caption").textContent = "";
  $("viewer").innerHTML = `<div class="viewer-empty">${message}</div>`;
  $("d-poses").innerHTML = `<span class="count">${list}</span>`;
}

function drawDetail() {
  const c = byId.get(state.selected);
  if (!c) return;
  $("d-id").textContent = c.id;
  $("d-depict").innerHTML = c.svg || "";
  $("d-props").innerHTML =
    `<span>Heavy atoms <b>${c.heavy}</b></span><span>Rot. bonds <b>${c.rotb}</b></span>` +
    `<span>MW <b>${c.mw}</b></span>` +
    (c.charge === undefined ? "" : `<span>Charge <b>${c.charge > 0 ? "+" : ""}${c.charge}</b></span>`) +
    `<span>Cofactors <b>${c.cofactors.length ? c.cofactors.join(", ") : "none"}</b></span>`;
  // IDs are <PDB code>_<ligand CCD code>; both have pages at the RCSB.
  const [pdb, ...rest] = c.id.split("_");
  const ccd = rest.join("_");
  const link = (href, text, title) => `<a href="${href}" target="_blank" rel="noopener" title="${title}">${text} ↗</a>`;
  $("d-protein").innerHTML =
    (c.residues ? `<span>Protein <b>${c.chains} chain${c.chains === 1 ? "" : "s"}, ${c.residues} residues</b></span>` : "") +
    `<span>${link(`https://www.rcsb.org/structure/${pdb}`, `PDB ${pdb}`, "Structure entry at the RCSB PDB")}</span>` +
    `<span>${link(`https://www.rcsb.org/ligand/${ccd}`, `Ligand ${ccd}`, "Ligand page at the RCSB PDB")}</span>`;
  $("d-smiles-row").hidden = !c.smiles;
  $("d-smiles").textContent = c.smiles || "";
  $("d-smiles").title = c.smiles || "";
  $("d-arms").innerHTML = ARMS.map((arm) => {
    const o = c.arms[arm.key];
    return (
      `<div class="arm-cell${arm.key === state.arm ? " sel" : ""}"><div class="name">${arm.label}${o.real ? " ◆" : ""}</div>` +
      `<div class="val"${o.error ? ` title="${o.error.replace(/"/g, "&quot;")}"` : ""}>${glyph(L.classify(o))}<span>${L.fmtA(o.rmsd)} Å</span>` +
      `<span class="best" title="Best RMSD among all returned poses">best ${L.fmtA(o.best)}</span></div></div>`
    );
  }).join("");

  if (!poseCache.has(c.id) && DATA.poseDir) {
    requestPoses(c.id);
    showNoPoses("Loading poses…", "Loading…");
    return;
  }
  const real = poseCache.get(c.id);
  const poses = armPoses(real);
  if (!real || !poses.length) {
    const o = c.arms[state.arm];
    const message = o.error
      ? `${armLabel(state.arm)} failed for this complex: ${o.error}`
      : DATA.demo
        ? "Docked poses for this complex will appear here after the full benchmark run.<br>Complexes marked ◆ have real poses now."
        : real === null
          ? "The pose file for this complex was not found. Rebuild the data with <code>scripts/explorer.sh</code>."
          : `No poses were saved for ${armLabel(state.arm)}.`;
    showNoPoses(message, "No poses.");
    return;
  }
  $("viewer-bar").hidden = false;
  $("poses-head").textContent = `${armLabel(Array.isArray(real.poses) ? BASE : state.arm)} poses · score · RMSD (Å)`;
  $("d-poses").innerHTML = poses
    .map((p) => {
      const ok = p.rmsd < 2 ? ` <span class="ok">✓</span>` : "";
      return (
        `<button type="button" data-rank="${p.rank}" aria-pressed="${p.rank === state.pose}" aria-label="Pose ${p.rank}, score ${p.score}, RMSD ${L.fmtA(p.rmsd)} Å">` +
        `<span class="rank">#${p.rank}</span><span>${p.score.toFixed(1)}</span><span>${p.rmsd.toFixed(2)}${ok}</span></button>`
      );
    })
    .join("");
  if (state.tab === "complexes") drawViewer(real, poses);
}

// Protein chains: muted hues, none blue, so the docked pose (blue carbons) stands out.
const CHAIN_COLOURS = ["#d1495b", "#b5b44a", "#3f9c8f", "#9a6cc8", "#e0883d", "#c76b98"];

/** Colour per chain ID, in order of first appearance in the PDB text. */
function chainColours(pdb) {
  const ids = [...new Set(pdb.split("\n").filter((l) => l.startsWith("ATOM")).map((l) => l[21]))];
  return Object.fromEntries(ids.map((id, i) => [id, CHAIN_COLOURS[i % CHAIN_COLOURS.length]]));
}

// The protein layers (cartoon, surface, pocket) are rebuilt only when the complex or a layer
// toggle changes; switching poses or arms swaps just the ligands, so the surface is computed once.
let scene = { key: null, ligands: [] };

function drawViewer(real, poses) {
  const el = $("viewer");
  if (typeof window.$3Dmol === "undefined") {
    el.innerHTML = `<div class="viewer-empty">The 3D viewer library did not load. Check your connection and reload.</div>`;
    return;
  }
  if (!viewer) {
    el.innerHTML = "";
    viewer = window.$3Dmol.createViewer(el, { backgroundColor: css("--viewer-bg"), antialias: true });
    scene = { key: null, ligands: [] };
  }
  const withProtein = state.showProtein && Boolean(real.protein);
  const key = [state.selected, withProtein, state.showSurface, state.showPocket, state.closeUp, css("--viewer-bg")].join("|");
  const rebuild = scene.key !== key;
  if (rebuild) {
    viewer.setBackgroundColor(css("--viewer-bg"));
    viewer.removeAllSurfaces();
    viewer.clear();
    scene = { key, ligands: [], protein: null };
    if (withProtein) {
      const map = chainColours(real.protein);
      scene.protein = viewer.addModel(real.protein, "pdb");
      // Opacity below 1 puts the cartoon in 3Dmol's transparent pass, where a surface hides it.
      // So: opaque cartoon inside the surface for the whole protein; in close-up, where the
      // cartoon sits in front of the ligand, a see-through cartoon and no surface.
      viewer.setStyle({ model: scene.protein }, { cartoon: { colorscheme: { prop: "chain", map }, opacity: state.closeUp ? 0.5 : 1 } });
      if (state.showSurface && !state.closeUp) {
        // One surface per chain: 3Dmol colours a surface from `color`, not from a chain map.
        for (const [chain, color] of Object.entries(map)) {
          const done = viewer.addSurface(window.$3Dmol.SurfaceType.MS, { opacity: 0.35, color }, { model: scene.protein, chain });
          if (done && typeof done.then === "function") done.then(() => viewer.render());
        }
      }
    }
    if (state.showPocket) {
      const pocket = viewer.addModel(real.pocket, "pdb");
      // Residues: thin and faded, so the two ligands stay in front. Cofactors: violet, whole.
      viewer.setStyle({ model: pocket }, { stick: { radius: 0.09, opacity: 0.75, colorscheme: { prop: "elem", map: elementColours(RESIDUE) } } });
      if (withProtein) viewer.setStyle({ model: pocket, atom: ["N", "C", "O"], hetflag: false }, {}); // backbone is in the cartoon
      viewer.setStyle({ model: pocket, hetflag: true }, { stick: { radius: 0.16, colorscheme: { prop: "elem", map: elementColours(COFACTOR) } } });
    }
  } else {
    scene.ligands.forEach((m) => viewer.removeModel(m));
  }
  const crystal = viewer.addModel(real.crystal, "sdf");
  viewer.setStyle({ model: crystal }, state.showCrystal ? { stick: { radius: 0.36, color: CRYSTAL, opacity: 0.7 } } : {});
  const pose = poses.find((p) => p.rank === state.pose) || poses[0];
  const docked = viewer.addModel(pose.block, "sdf");
  viewer.setStyle({ model: docked }, { stick: { radius: 0.18, colorscheme: { prop: "elem", map: elementColours(css("--metric-1")) } } });
  scene.ligands = [crystal, docked];
  viewer.resize();
  if (rebuild) {
    // Whole protein by default, as in a structure figure. Close-up keeps the protein's depth
    // range (so the cartoon is not clipped away) and moves in on the ligand.
    if (scene.protein) {
      viewer.zoomTo({ model: scene.protein });
      if (!state.closeUp) viewer.zoom(1.3); // zoomTo leaves a wide margin
      if (state.closeUp) {
        viewer.center({ model: crystal });
        viewer.zoom(3.5);
      }
    } else {
      viewer.zoomTo({ model: crystal });
      viewer.zoom(0.7);
    }
  }
  viewer.render();
  const ok = pose.rmsd < 2;
  $("viewer-caption").innerHTML =
    `<b>Blue</b>: docked pose #${pose.rank} (${armLabel(Array.isArray(real.poses) ? BASE : state.arm)}). <b>Amber</b>: crystal pose, the right answer. ` +
    `They are ${L.fmtA(pose.rmsd)} Å apart (RMSD), ${ok ? "under 2 Å: <span class=\"delta up\">a success</span>" : "over 2 Å: <span class=\"delta down\">a miss</span>"}.`;
}

// The reference pose: one solid colour, drawn thicker than the docked pose so that a close
// match shows as an amber halo around the blue pose rather than hiding behind it.
const CRYSTAL = "#f08c00";
const RESIDUE = "#a7afbb";
const COFACTOR = "#8e6cc8";

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

/** Copy text to the clipboard, confirming on the button; shows the text if copying is refused. */
async function copyText(button, text, label) {
  try {
    await navigator.clipboard.writeText(text);
    button.textContent = "Copied";
  } catch {
    button.textContent = text;
  }
  setTimeout(() => { button.textContent = label; }, 2500);
}

/** Show the selected ligand's 2D structure large, with its full SMILES. */
function openDepiction() {
  const c = byId.get(state.selected);
  if (!c || !c.svg) return;
  $("depict-title").textContent = c.id;
  $("depict-large").innerHTML = c.svg;
  $("depict-smiles").textContent = c.smiles || "";
  const dialog = $("depict-dialog");
  if (typeof dialog.showModal === "function") dialog.showModal();
  else dialog.setAttribute("open", ""); // older browsers and jsdom
}

function copyLink() {
  const base = DATA.shareUrl || `${location.origin}${location.pathname}`;
  const engine = ENGINES ? `?engine=${state.engine}` : "";
  return copyText($("copy-link"), `${base}${engine}#${state.selected}`, "Copy link");
}

function fillArmSelect() {
  $("arm-select").innerHTML = ARMS.filter((a) => a.key !== BASE)
    .map((a) => `<option value="${a.key}">${a.label}: ${a.detail}</option>`)
    .join("");
  $("arm-select").value = state.arm;
}

function showMeta() {
  $("meta").textContent = `${DATA.dataset} · ${DATA.complexes.length} complexes · ${DATA.engine}`;
}

function setEngine(key) {
  state.engine = useEngine(key);
  if (!ARMS.some((a) => a.key === state.arm)) state.arm = defaultArm(ARMS);
  state.page = 1;
  fillArmSelect();
  showMeta();
  try {
    history.replaceState(null, "", `?engine=${state.engine}${location.hash}`);
  } catch {
    /* sandboxed frames may refuse */
  }
  renderAll();
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
  showMeta();
  $("engine-field").hidden = !ENGINES;
  if (ENGINES) {
    $("engine-select").innerHTML = ENGINES.map((e) => `<option value="${e.key}">${e.label}</option>`).join("");
    $("engine-select").value = state.engine;
    $("engine-select").addEventListener("change", (e) => setEngine(e.target.value));
  }
  $("badge").hidden = !DATA.demo;
  $("real-key").hidden = !DATA.demo;
  $("key").innerHTML = CLASSES.map((k) => `<span>${glyph(k)} ${CLASS_LABEL[k]}</span>`).join("");
  $("footer").textContent =
    "Outcome classes follow the top-ranked pose: valid success (RMSD < 2 Å and PB-valid); correct but not PB-valid; " +
    `scoring miss (a pose under 2 Å exists below rank 1); wrong site (centroid more than ${SITE_SHIFT} Å away, a heuristic); ` +
    "wrong binding mode (right site, wrong pose); failed (preparation or docking raised an error, counted as a miss).";
  fillArmSelect();

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
  $("t-protein").addEventListener("change", (e) => { state.showProtein = e.target.checked; drawDetail(); });
  $("t-surface").addEventListener("change", (e) => { state.showSurface = e.target.checked; drawDetail(); });
  $("t-closeup").addEventListener("change", (e) => {
    state.closeUp = e.target.checked;
    $("t-surface").disabled = state.closeUp; // no surface in close-up: it would hide the ligand
    drawDetail();
  });
  $("copy-link").addEventListener("click", copyLink);
  $("d-depict").addEventListener("click", openDepiction);
  $("depict-close").addEventListener("click", () => $("depict-dialog").close());
  $("depict-dialog").addEventListener("click", (e) => {
    if (e.target === $("depict-dialog")) $("depict-dialog").close(); // click on the backdrop
  });
  $("copy-smiles").addEventListener("click", () => copyText($("copy-smiles"), byId.get(state.selected).smiles, "Copy"));

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
