# Changelog

All notable changes to **LINwalker** will be documented in this file.

The format is based on *Keep a Changelog*, and this project follows *Semantic Versioning*.

## 1.1.0
### Added
- **place:** reference-anchored placement of uploaded/unlabelled cgMLST profiles against official LIN-coded PubMLST genomes.
- Nearest-reference ST and clonal-complex context, plus exact official cgST reporting only for complete zero-distance reference-profile matches.
- Robust parsing of BIGSdb Genome Comparator uploaded-genome labels (`uN (...)`) back to stable query IDs.
- Direct parsing of PubMLST Genome Comparator Excel `all` worksheets.
- BIGSdb-style normalisation of allele differences when loci are missing in either genome.
- Conservative LIN-prefix placement using configurable support/proportion thresholds.
- Conservative `Cjc_cgc2_200/100/50/25/10/5` placement where labelled reference groups are available.
- Pairwise, threshold-by-threshold, cgc2 and reference-label audit tables.
- Tests for current 18-level Campylobacter cgMLST-v2 LINcodes and missing-locus normalisation.

### Changed
- CLI LIN depth can now be inferred from the supplied LINcodes instead of assuming 17 levels.
- Current Campylobacter cgMLST-v2 worked examples use the 18-component LIN hierarchy.

### Interpretation
- `place` does **not** create official PubMLST cgSTs or LINcodes. It reports nearest official anchors and conservative supported placement within the existing PubMLST hierarchy.

## 1.0.16
### Fixed
- **stcc:** Coerce purity columns to numeric before plotting, preventing matplotlib crashes on `pandas.NA`/`NAType`.
- **outbreak:** Drop missing/invalid LIN codes before building LIN prefixes, avoiding a dominant `nan_nan_...` cluster.

## 1.0.17
### Fixed
- **stcc:** More robust plotting (drop-NA per series) and clearer behaviour when ST/CC columns are missing (annotated figure instead of "empty" plot).
- **stcc:** Auto-detect common PubMLST ST / clonal complex column headers when the default names aren't present.
- **outbreak:** Exclude prefixes containing literal missing tokens (e.g. `nan_nan_...`) even when they appear as strings.

## 1.0.18
### Added
- **diversify:** Rarefied diversification curves (sample-size normalisation) for reservoir-only and all-sources outputs.
### Changed
- **diversify:** X-axis tick labels simplified to show only LIN thresholds 1, 5, 10, and 15.

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

