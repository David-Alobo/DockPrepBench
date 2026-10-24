// @vitest-environment jsdom
/**
 * Integration tests: load the real index.html and app.js into jsdom with fixture data and
 * drive the page as a user would. The 3D viewer library is absent here, which exercises the
 * page's fallback.
 */
import { readFileSync } from "node:fs";
import { createRequire } from "node:module";
import { dirname, join } from "node:path";
import { fileURLToPath } from "node:url";
import { beforeEach, describe, expect, it } from "vitest";

import { makeData } from "./fixture.js";

const here = dirname(fileURLToPath(import.meta.url));
const require = createRequire(import.meta.url);
const html = readFileSync(join(here, "..", "index.html"), "utf8");
const appSource = readFileSync(join(here, "..", "app.js"), "utf8");
const body = html.match(/<body>([\s\S]*)<\/body>/)[1].replace(/<script[\s\S]*?<\/script>/g, "");

const $ = (id) => document.getElementById(id);
const text = (id) => $(id).textContent.trim();
const visibleRows = () => [...document.querySelectorAll("#rows tr[data-id]")].map((tr) => tr.dataset.id);
const click = (el) => el.dispatchEvent(new MouseEvent("click", { bubbles: true }));

/** Start the page fresh, optionally opening a link such as "#C23_LIG". */
function boot(hash = "") {
  document.body.innerHTML = body;
  window.history.replaceState(null, "", `/${hash}`);
  window.EXPLORER_DATA = makeData();
  window.ExplorerLib = require("../lib.js");
  delete window.$3Dmol;
  window.matchMedia = () => ({ matches: false, addEventListener() {} });
  new Function(appSource)();
}

describe("overview", () => {
  beforeEach(() => boot());

  it("opens on the overview tab", () => {
    expect($("panel-overview").hidden).toBe(false);
    expect($("panel-complexes").hidden).toBe(true);
    expect($("tab-overview").getAttribute("aria-selected")).toBe("true");
  });

  it("shows the comparison tiles for the default arm", () => {
    expect(text("t-base")).toBe("40.0%");
    expect(text("t-arm")).toBe("56.0%");
    expect(text("t-arm-ci")).toMatch(/^\+16\.0 pts vs baseline/);
    expect(text("t-changed")).toBe("6");
    expect(text("t-changed-s")).toBe("↑ 5 gained · ↓ 1 lost");
  });

  it("updates every tile when the success metric changes", () => {
    $("metric-select").value = "valid";
    $("metric-select").dispatchEvent(new Event("change"));
    expect(text("t-base")).toBe("36.0%");
    expect(text("t-arm")).toBe("52.0%");
  });

  it("switches arms from the dropdown", () => {
    $("arm-select").value = "ligand_dimorphite";
    $("arm-select").dispatchEvent(new Event("change"));
    expect(text("t-arm")).toBe("40.0%");
    expect(text("t-changed")).toBe("0");
    expect(text("t-p")).toBe("p = 1.0");
  });

  it("draws one chart row per arm, baseline first", () => {
    const labels = [...document.querySelectorAll("#chart .arm-label")].map((t) => t.textContent);
    expect(labels).toEqual(["Baseline (ref.)", "Box 2 Å", "Protonated ligand"]);
  });

  it("fills the transition matrix with counts that sum to the dataset", () => {
    const cells = [...document.querySelectorAll("#matrix td")];
    expect(cells).toHaveLength(36); // six classes, failed included
    expect(cells.reduce((s, td) => s + Number(td.textContent), 0)).toBe(25);
    expect(document.querySelector('#matrix td[data-cell="scoring|valid"]').textContent).toBe("5");
  });
});

describe("complexes tab", () => {
  beforeEach(() => {
    boot();
    click($("tab-complexes"));
  });

  it("shows ten complexes per page with page controls", () => {
    expect(visibleRows()).toHaveLength(10);
    expect(text("count")).toBe("1–10 of 25");
    expect(text("page-label")).toBe("1 / 3");
    expect($("page-prev").disabled).toBe(true);

    click($("page-next"));
    expect(visibleRows()[0]).toBe("C10_LIG");
    expect(text("count")).toBe("11–20 of 25");

    click($("page-last"));
    expect(visibleRows()).toEqual(["C20_LIG", "C21_LIG", "C22_LIG", "C23_LIG", "C24_LIG"]);
    expect(document.querySelectorAll("#rows tr")).toHaveLength(10); // filler keeps the height
    expect($("page-next").disabled).toBe(true);
  });

  it("returns to page 1 when a search changes, even if later pages still exist", () => {
    click($("page-last"));
    $("search").value = "_lig"; // still matches all 25, so page 3 would remain valid
    $("search").dispatchEvent(new Event("input"));
    expect(text("page-label")).toBe("1 / 3");
    expect(visibleRows()[0]).toBe("C00_LIG");
  });

  it("returns to page 1 when the sort changes", () => {
    click($("page-next"));
    click(document.querySelector('th button[data-sort="rotb"]'));
    expect(text("page-label")).toBe("1 / 3");
  });

  it("reports filtered counts", () => {
    $("search").value = "c1";
    $("search").dispatchEvent(new Event("input"));
    expect(visibleRows()).toHaveLength(10);
    expect(text("count")).toBe("1–10 of 10 (filtered from 25)");
  });

  it("says so when nothing matches", () => {
    $("search").value = "nothing";
    $("search").dispatchEvent(new Event("input"));
    expect(visibleRows()).toHaveLength(0);
    expect(document.querySelector(".empty-row").textContent).toMatch(/No complexes match/);
  });

  it("sorts by a column and toggles direction", () => {
    const heavy = document.querySelector('th button[data-sort="heavy"]');
    click(heavy);
    expect(visibleRows()[0]).toBe("C00_LIG");
    click(heavy);
    expect(visibleRows()[0]).toBe("C24_LIG");
    expect(heavy.dataset.active).toBe("-1");
  });

  it("shows the selected complex and falls back when the 3D library is missing", () => {
    click(document.querySelector('#rows tr[data-id="C03_LIG"]'));
    expect(text("d-id")).toBe("C03_LIG");
    expect(document.querySelector('#rows tr[aria-selected="true"]').dataset.id).toBe("C03_LIG");
    expect($("d-arms").querySelectorAll(".arm-cell")).toHaveLength(3);
    expect($("viewer").textContent).toMatch(/will appear here/);
  });
});

describe("built data: failures and poses loaded per complex", () => {
  /** Boot with data as build_data.py writes it: no demo, a failed complex, a pose directory. */
  function bootBuilt() {
    document.head.querySelectorAll('script[src^="poses/"]').forEach((s) => s.remove()); // from earlier tests
    boot();
    const data = makeData();
    data.demo = false;
    data.poseDir = "poses/test";
    data.complexes[24].arms.box_padding_2 = { error: "RuntimeError: receptor failed" };
    Object.assign(data.complexes[3], { charge: -1, chains: 2, residues: 300 });
    window.EXPLORER_DATA = data;
    document.body.innerHTML = body;
    new Function(appSource)();
    click($("tab-complexes"));
  }
  const payload = {
    crystal: "",
    pocket: "",
    poses: {
      baseline: [{ rank: 1, score: -7.1, rmsd: 0.8, block: "" }],
      box_padding_2: [
        { rank: 1, score: -8.2, rmsd: 3.4, block: "" },
        { rank: 2, score: -8.0, rmsd: 1.1, block: "" },
      ],
    },
  };

  beforeEach(bootBuilt);

  it("shows the ligand SMILES, charge, protein size and RCSB links", () => {
    window.EXPLORER_DATA.complexes[3].smiles = "CC(=O)O";
    click(document.querySelector('#rows tr[data-id="C03_LIG"]'));
    expect(text("d-smiles")).toBe("CC(=O)O");
    expect($("d-smiles-row").hidden).toBe(false);
    expect(text("d-props")).toMatch(/Charge -1/);
    expect(text("d-protein")).toMatch(/2 chains, 300 residues/);
    const links = [...$("d-protein").querySelectorAll("a")].map((a) => a.getAttribute("href"));
    expect(links).toEqual(["https://www.rcsb.org/structure/C03", "https://www.rcsb.org/ligand/LIG"]);
  });

  it("hides the prototype badge and the ◆ key", () => {
    expect($("badge").hidden).toBe(true);
    expect($("real-key").hidden).toBe(true);
  });

  it("shows a failed complex as failed, without a distance", () => {
    click($("page-last"));
    click(document.querySelector('#rows tr[data-id="C24_LIG"]'));
    const cell = [...$("d-arms").querySelectorAll(".arm-cell")][1];
    expect(cell.textContent).toMatch(/– Å/);
    expect(document.querySelector('#matrix td[data-cell="mode|failed"]').textContent).toBe("1");
  });

  it("requests a complex's pose file once, then shows the selected arm's poses", () => {
    const scripts = () => [...document.head.querySelectorAll('script[src^="poses/test/"]')];
    click(document.querySelector('#rows tr[data-id="C03_LIG"]'));
    expect(scripts().map((s) => s.getAttribute("src"))).toContain("poses/test/C03_LIG.js");
    expect($("viewer").textContent).toMatch(/Loading/);

    window.ExplorerPoses("C03_LIG", payload);
    expect($("d-poses").querySelectorAll("button[data-rank]")).toHaveLength(2); // box_padding_2
    expect(text("poses-head")).toMatch(/^Box 2 Å poses/);
    click(document.querySelector('#rows tr[data-id="C03_LIG"]'));
    expect(scripts().filter((s) => s.getAttribute("src") === "poses/test/C03_LIG.js")).toHaveLength(1);
  });

  it("explains a missing pose file", () => {
    click(document.querySelector('#rows tr[data-id="C04_LIG"]'));
    window.ExplorerPoses("C04_LIG", null);
    expect($("viewer").textContent).toMatch(/not found/);
  });
});

describe("several engines", () => {
  /** Two engines: A is the fixture; B has only box_padding_2, identical to its baseline. */
  function bootEngines(search = "") {
    document.head.querySelectorAll('script[src^="poses/"]').forEach((s) => s.remove());
    document.body.innerHTML = body;
    window.history.replaceState(null, "", `/${search}`);
    const data = makeData();
    const armsB = data.arms.slice(0, 2);
    data.engines = [
      { key: "a", label: "Engine A", arms: data.arms },
      { key: "b", label: "Engine B", arms: armsB },
    ];
    for (const c of data.complexes) {
      c.byEngine = { a: c.arms, b: { baseline: c.arms.baseline, box_padding_2: c.arms.baseline } };
      delete c.arms;
    }
    data.demo = false;
    data.poseDir = "poses";
    delete data.arms;
    window.EXPLORER_DATA = data;
    window.ExplorerLib = require("../lib.js");
    delete window.$3Dmol;
    window.matchMedia = () => ({ matches: false, addEventListener() {} });
    new Function(appSource)();
  }
  const pick = (value) => {
    $("engine-select").value = value;
    $("engine-select").dispatchEvent(new Event("change"));
  };

  it("starts on the first engine and lists every engine", () => {
    bootEngines();
    expect($("engine-field").hidden).toBe(false);
    expect([...$("engine-select").options].map((o) => o.textContent)).toEqual(["Engine A", "Engine B"]);
    expect(text("meta")).toMatch(/Engine A$/);
    expect(text("t-changed")).toBe("6");
  });

  it("switches outcomes, arms and the header with the engine", () => {
    bootEngines();
    pick("b");
    expect(text("meta")).toMatch(/Engine B$/);
    expect(text("t-changed")).toBe("0");
    expect([...$("arm-select").options].map((o) => o.value)).toEqual(["box_padding_2"]);
    expect(window.location.search).toBe("?engine=b");
  });

  it("falls back to a default arm when the new engine lacks the selected one", () => {
    bootEngines();
    $("arm-select").value = "ligand_dimorphite";
    $("arm-select").dispatchEvent(new Event("change"));
    pick("b");
    expect($("arm-select").value).toBe("box_padding_2");
  });

  it("opens the engine named in the link", () => {
    bootEngines("?engine=b");
    expect($("engine-select").value).toBe("b");
    expect(text("meta")).toMatch(/Engine B$/);
  });

  it("shows the selected engine's poses", () => {
    bootEngines();
    click($("tab-complexes"));
    click(document.querySelector('#rows tr[data-id="C03_LIG"]'));
    const pose = (rank, rmsd) => ({ rank, score: -7, rmsd, block: "" });
    window.ExplorerPoses("C03_LIG", {
      crystal: "", pocket: "",
      poses: { a: { box_padding_2: [pose(1, 1), pose(2, 3)] }, b: { box_padding_2: [pose(1, 5)] } },
    });
    expect($("d-poses").querySelectorAll("button[data-rank]")).toHaveLength(2);
    pick("b");
    expect($("d-poses").querySelectorAll("button[data-rank]")).toHaveLength(1);
  });
});

describe("links between views", () => {
  it("opens the complexes behind a matrix cell, and clears that filter", () => {
    boot();
    click(document.querySelector('#matrix td[data-cell="scoring|valid"]'));
    expect($("panel-complexes").hidden).toBe(false);
    expect(visibleRows()).toEqual(["C10_LIG", "C11_LIG", "C12_LIG", "C13_LIG", "C14_LIG"]);
    expect($("active-filter").hidden).toBe(false);
    expect(text("d-id")).toBe("C10_LIG");

    click($("clear-cell"));
    expect($("active-filter").hidden).toBe(true);
    expect(text("count")).toBe("1–10 of 25");
  });

  it("opens a linked complex on its own page", () => {
    boot("#c23_lig");
    expect($("panel-complexes").hidden).toBe(false);
    expect(text("d-id")).toBe("C23_LIG");
    expect(text("page-label")).toBe("3 / 3");
  });
});
