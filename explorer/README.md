# Results explorer

A web page for browsing benchmark results. It shows success by experiment arm, how
outcomes moved between the baseline and an arm, and each complex with a 3D view of the
crystal ligand and the docked poses in its binding pocket. It is plain HTML, CSS and
JavaScript with no build step. [3Dmol.js](https://3dmol.org) and the IBM Plex fonts load
from CDNs.

## Start it

```bash
scripts/explorer.sh                       # build the data for Vina and serve the page
scripts/explorer.sh --engine gnina --open # another engine; open the browser too
```

Then open <http://localhost:8000/>. The script needs [uv](https://docs.astral.sh/uv/) and
works on Linux, macOS and Git Bash on Windows. Options:

| Option | Default | Meaning |
|---|---|---|
| `--engine` | `vina` | Engine whose runs to show (`vina`, `smina`, `gnina`, `rdock`) |
| `--results DIR` | `results/server/results` if present, else `results` | Root holding `<engine>/<arm>/` runs |
| `--port N` | `8000` | Port to serve on (bound to localhost only) |
| `--skip-build` | off | Serve the data built last time |
| `--open` | off | Open the default browser |
| `--share-url URL` | page address | Base URL used by "Copy link", e.g. a hosted copy |

Every finished arm of the engine is shown; seed replicates are left out. A complex whose
preparation or docking failed is shown as *failed* and counted as a miss, as in the paper.

## Build the data only

```bash
uv run python explorer/build_data.py --engine vina
```

This writes `data.js` (outcomes and ligand properties, about 2 MB) and
`poses/<engine>/<id>.js` (crystal ligand, binding pocket from the dataset's protein file,
and every arm's docked poses; about 150 KB each). The page loads a complex's pose file
only when the complex is opened. Both are generated and ignored by git. The PoseBusters
set is downloaded to `data/raw/` on first use.

Opening `index.html` straight from the file system also works once the data is built.

Links to a complex are the page address followed by `#<complex id>`, e.g.
`http://localhost:8000/#5S8I_2LY`. They open that complex on the Complexes tab.

## Tests

Unit tests cover the logic in `lib.js` (outcome classes including failures, statistics
checked against the Python implementation, filtering, paging). Integration tests load
`index.html` and `app.js` into jsdom with a fixture dataset and drive the page: tabs,
tiles, pagination, sorting, the matrix-to-table link, complex links, failed complexes and
on-demand pose loading. They run in CI.

```bash
cd explorer
npm ci
npm test            # or: npm run test:watch
```

## Files

| File | Role |
|---|---|
| `index.html` | Page structure and styles (light and dark themes) |
| `lib.js` | Pure logic: outcome classes, statistics (Wilson intervals, exact McNemar), filtering, paging |
| `app.js` | Page wiring: chart, transition matrix, table, detail panel, 3D viewer |
| `test/` | Vitest unit and jsdom integration tests |
| `build_data.py` | Writes `data.js` and `poses/` from the dataset and a results directory |
| `data.js` | Generated data (`window.EXPLORER_DATA`), not committed |
| `poses/` | Generated per-complex pose files, loaded on demand, not committed |
