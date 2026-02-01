# Changelog

## v1.0.14

- **Fix**: diversification plot x-axis tick labels are no longer squashed; major tick labels shown at **1, 5, 10, 15** with minor ticks for each level.
- **Fix**: CLI/plotting alignment across modules (introgress, stcc, tree, outbreak) using a consistent `outdir/plots|tables|logs` structure.
- **Add**: stable output structure created automatically (`plots/`, `tables/`, `logs/`, plus `derived/` for prep outputs).
- **Add**: `pyproject.toml` so `pip install -e .` works.
- **Add**: `CITATION.cff` and CI smoke tests.

## v1.0.13 and earlier

- Development iterations during initial packaging and CLI expansion.
