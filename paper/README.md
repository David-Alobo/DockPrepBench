# Manuscript

LaTeX source for the DockPrepBench paper, intended first as a ChemRxiv preprint. It uses a
journal-neutral `article` layout with numbered citations, so it can move to a journal
template (e.g. `achemso` for ACS journals) without rewriting.

## Build

The PDF is built with `latexmk`, which runs pdflatex and BibTeX as often as needed. It
writes everything to `paper/build/` (git-ignored).

```bash
# in WSL Ubuntu (TeX Live installed), from the repository root
cd paper && latexmk main.tex        # → paper/build/main.pdf
latexmk -c                          # remove intermediate files
```

CI builds the PDF whenever `paper/` or `docs/figures/` changes (workflow `Paper`) and fails on LaTeX errors or undefined citations or references. Download the current PDF from the workflow run's `manuscript` artifact on GitHub (Actions → Paper → latest run).

On Overleaf, upload `paper/` plus the figures it uses from `docs/figures/` (adjust
`\graphicspath` in `main.tex` if the folder layout changes).

## Layout

| Path | Content |
|---|---|
| `main.tex` | Preamble, title, section order, macros |
| `sections/` | One file per section; Methods drafted, others outlined |
| `references.bib` | Bibliography; every journal entry checked against Crossref |
| `latexmkrc` | Build settings |

Figures come from `docs/figures/` (study design) and `figures/` (made by
`dockprepbench figures`), so the paper always uses the current versions.

## Conventions

- `\result{...}` (orange) marks a number or statement that depends on runs not yet done;
  `\todo{...}` (red) marks writing still to do. Both must be gone before submission:
  `grep -rn "\\\\result{\|\\\\todo{" paper/sections`.
- Methods must match the code and `docs/methods.md`: change all three in the same pull request.
- Report results following the `benchmark-reporting` conventions: denominators, both
  success metrics, 95% Wilson intervals, effect size before the p-value, Holm-adjusted
  p-values across arms.

## References

Add a reference only after checking its DOI, author list, journal, volume, pages and year
against Crossref (`https://api.crossref.org/works/<doi>`) or the publisher's page. Record
candidates not yet checked here rather than in `references.bib`.

To verify before citing:

- Prior studies of protonation, tautomer and water effects on docking, for the
  Introduction (to be found).
