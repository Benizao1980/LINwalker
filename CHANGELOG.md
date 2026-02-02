# Changelog

All notable changes to **LINwalker** will be documented in this file.

The format is based on *Keep a Changelog*, and this project follows *Semantic Versioning*.

## 1.0.16
### Fixed
- **stcc:** Coerce purity columns to numeric before plotting, preventing matplotlib crashes on `pandas.NA`/`NAType`.
- **outbreak:** Drop missing/invalid LIN codes before building LIN prefixes, avoiding a dominant `nan_nan_...` cluster.

## 1.0.17
### Fixed
- **stcc:** More robust plotting (drop-NA per series) and clearer behaviour when ST/CC columns are missing (annotated figure instead of "empty" plot).
- **stcc:** Auto-detect common PubMLST ST / clonal complex column headers when the default names aren't present.
- **outbreak:** Exclude prefixes containing literal missing tokens (e.g. `nan_nan_...`) even when they appear as strings.

## 1.0.15
### Added
- Restored full CLI command set: `prep`, `diversify`, `introgress`, `stcc`, `tree`, `outbreak`.
- Package metadata for editable installs (`pyproject.toml`).

### Changed
- Standardised plot output handling to support `--formats png svg` across modules.

## 1.0.14
### Added
- CI smoke tests scaffold.
- Stable output structure conventions (`plots/`, `tables/`, `logs/`).

