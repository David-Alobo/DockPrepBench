# Results explorer (prototype)

A static web page for browsing benchmark results. It shows success by experiment arm, how
outcomes moved between the baseline and an arm, and each complex with a 3D view of
crystal vs docked poses. It is plain HTML, CSS and JavaScript with no build step or
server. [3Dmol.js](https://3dmol.org) and the IBM Plex fonts load from CDNs.

**Status: prototype.** Arm outcomes in `data.js` are simulated until the benchmark runs
finish. The complex IDs, ligand properties, cofactors, and the baseline poses of
complexes with a real run are real. The builder will switch to real results from
`results/<engine>/<arm>/`.

## Build the data

`data.js` is generated, and git ignores it:

```bash
uv run python explorer/build_data.py . explorer
```

An optional third argument sets the base URL used by "Copy link" (for example the hosted
page). Without it the page uses its own address.

## View it

Open `explorer/index.html` in a browser. It works from the file system, without a
server. An internet connection is needed for the 3D viewer library and the fonts.

Links to a complex are the page address followed by `#<complex id>`, e.g.
`index.html#5S8I_2LY`. They open that complex on the Complexes tab.

## Tests

Unit tests cover the logic in `lib.js` (outcome classes, statistics checked against the
Python implementation, filtering, paging). Integration tests load `index.html` and `app.js`
into jsdom with a fixture dataset and drive the page: tabs, tiles, pagination, sorting,
the matrix-to-table link and complex links. They run in CI.

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
| `build_data.py` | Writes `data.js` from the dataset, results and saved poses |
| `data.js` | Generated data (`window.EXPLORER_DATA`), not committed |
