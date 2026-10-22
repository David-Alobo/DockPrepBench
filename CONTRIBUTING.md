# Contributing

## Setup

```bash
uv sync
uv run pre-commit install
```

## Workflow

1. Create a branch from `main`.
2. Make your change with tests in `tests/` mirroring the module under `src/dockprepbench/`.
3. Run the checks:

   ```bash
   uv run ruff check .
   uv run ruff format .
   uv run pytest
   ```

4. Add an entry under `[Unreleased]` in `CHANGELOG.md`.
5. Open a pull request; CI must pass on Linux, Windows and macOS.

## Conventions

- Raise exceptions on invalid input rather than returning `None`; silent failures corrupt benchmark results.
- Large inputs and generated outputs belong in `data/` and `results/`, which are not committed.
