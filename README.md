# LINscope

**LINscope** is a lightweight Python module for exploring hierarchical population structure and interspecies introgression using LIN codes derived from cgMLST profiles.

LINscope is designed for pathogens with extensive recombination, where population structure and ancestry vary across genomic scales. Rather than treating LIN codes as flat identifiers, LINscope explicitly exploits their hierarchical structure to quantify how genomic similarity, host association, and species identity change with resolution.

The module is intentionally minimal, deterministic, and analysis-focused.

---

## Core concepts

LIN codes encode hierarchical genomic relatedness across ordered levels, from coarse lineage structure to fine-scale genomic similarity.

LINscope uses this hierarchy to ask three biologically motivated questions:

1. **At what genomic scales does population structure emerge?**
2. **Where does species identity break down due to introgression?**
3. **How are these patterns modulated by host ecology?**

---

## Key functionality

### 1. LIN diversification by scale

Quantifies how rapidly genomic diversity accumulates across LIN thresholds by counting the number of unique truncated LIN identifiers at each hierarchical level.

Typical use cases:
- identifying scale-dependent host-associated structure
- selecting biologically meaningful LIN thresholds for downstream analyses
- comparing population complexity between reservoirs

---

### 2. Mixed-species LIN detection

Identifies LIN prefixes that contain isolates from multiple species and quantifies the proportion of mixed-species LIN clusters across hierarchical thresholds.

This provides a scale-aware view of interspecies introgression that is robust to mosaic genomes and recombination.

---

### 3. LIN Species Discordance Depth (LSDD)

Defines a per-isolate metric describing the earliest LIN level at which the species composition of an isolate’s LIN cluster diverges from its assigned species.

Lower values indicate deeper introgression; higher values indicate shallow or absent introgression.

---

## Design principles

- **Hierarchy-aware**: LIN codes are treated as ordered, multi-resolution descriptors
- **Deterministic**: no stochastic models or classifiers
- **Minimal scope**: focused on interpretation, not prediction
- **Composable**: intended to integrate with existing pipelines (e.g. cgMLST, SourceRunnerML)

---

## Intended use

LINscope is not a standalone pipeline and does not perform genome assembly, typing, or source attribution directly.

It is intended to:
- support exploratory population genomic analyses
- inform threshold selection for LIN-based classification
- provide interpretable summaries of recombination and introgression
- generate publication-ready figures

---

## Status

LINscope is under active development.

Initial applications focus on *Campylobacter jejuni* and *C. coli*, but the approach is general and applicable to any organism with hierarchical LIN codes.

---

## License

MIT
