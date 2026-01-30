# LINwalker

**LINwalker** is a lightweight Python module for exploring hierarchical population structure and interspecies introgression using LIN codes derived from cgMLST profiles.

LINwalker treats LIN codes as ordered, multi-resolution descriptors rather than flat identifiers, enabling scale-aware analysis of host association, lineage structure, and recombination.

The module is intentionally minimal, deterministic, and focused on interpretation rather than prediction.

---

## Core ideas

LIN codes encode hierarchical genomic similarity across ordered levels, from coarse lineage structure to near-isolate resolution.

LINwalker exploits this hierarchy to address three biological questions:

1. At what genomic scales does population structure emerge?
2. Where does species identity break down due to introgression?
3. How does host ecology modulate these patterns?

---

## Key functionality

### 1. LIN diversification by scale

Counts the number of unique truncated LIN identifiers across hierarchical thresholds to quantify how genomic diversity accumulates with resolution.

Typical uses:
- identifying host-associated structure
- selecting biologically meaningful LIN thresholds
- comparing population complexity between reservoirs

---

### 2. Mixed-species LIN detection

Identifies LIN prefixes that contain isolates from multiple species and quantifies the fraction of mixed-species LIN clusters across thresholds.

This provides a scale-aware view of introgression that is robust to recombination and mosaic genomes.

---

### 3. LIN Species Discordance Depth (LSDD)

Defines a per-isolate metric describing the earliest LIN level at which the species composition of an isolate’s LIN cluster diverges from its assigned species.

Lower values indicate deeper introgression; higher values indicate shallow or absent introgression.

---

## Design principles

- Hierarchy-aware: LIN codes are treated as ordered descriptors
- Deterministic: no stochastic models or classifiers
- Minimal scope: focused on interpretation, not pipelines
- Composable: designed to integrate with existing workflows

---

## Intended use

LINwalker is not a standalone typing or source attribution pipeline.

It is intended to:
- support exploratory population genomic analyses
- inform threshold selection for LIN-based classification
- quantify interspecies introgression
- generate publication-ready figures

---

## Status

LINwalker is under active development.

Initial applications focus on *Campylobacter jejuni* and *C. coli*, but the approach is general and applicable to any organism with hierarchical LIN codes.

---

## License

MIT
