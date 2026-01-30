# LINwalker

**LINwalker** is a lightweight Python toolkit for exploring hierarchical population structure and interspecies introgression using LIN codes (e.g., from PubMLST cgMLST schemes).

LINwalker treats LIN codes as **ordered, multi-resolution descriptors** rather than flat identifiers, enabling scale-aware analysis of host association, lineage structure, and recombination/introgression.

## What it does

LINwalker provides three core analyses:

1. **Diversification by scale** — how the number of unique LIN IDs grows as you increase LIN resolution.
2. **Mixed-species LINs** — how often LIN prefixes contain more than one species (scale-aware introgression signal).
3. **LSDD (LIN Species Discordance Depth)** — per-isolate metric for the earliest LIN level at which its LIN cluster’s majority species differs from its assigned species.

## Installation

### Option A: pip (recommended for most users)

```bash
pip install -r requirements.txt
```

### Option B: conda

```bash
conda create -n linwalker python=3.11 -y
conda activate linwalker
pip install -r requirements.txt
```

## Quick start

Minimal input is a tab-delimited table with at least:

- `LINcode (...)` column (underscore-separated levels)
- `source` (host/reservoir; optional for introgression-only workflows)
- `species` (required for mixed-species/LSDD analyses)

## Preparing PubMLST exports

PubMLST cgMLST exports often contain thousands of columns, many of which are not required for LIN-based analyses.

LINwalker provides a lightweight helper to convert raw PubMLST exports into analysis-ready tables:

```python
from linwalker.prep import prepare_pubmlst_export

lin_df, meta_df, cg_df = prepare_pubmlst_export(
    "PATHSAFE_pubmlst_export.tsv.gz"
)
```

Example usage:

```python
import pandas as pd
from linwalker.diversification import lin_diversification
from linwalker.introgression import mixed_species_summary, lsdd
from linwalker.plotting import plot_diversification

df = pd.read_csv("PATHSAFE_LINwalker_min.tsv", sep="\t")

# Diversification curves (unique LINs vs threshold) by source
div = lin_diversification(
    df,
    lin_col="LINcode (C. jejuni / C. coli cgMLST v2)",
    group_col="source"
)
plot_diversification(div, title="LIN diversification by host reservoir")

# Mixed-species LINs across thresholds (requires df['species'])
mix = mixed_species_summary(
    df,
    lin_col="LINcode (C. jejuni / C. coli cgMLST v2)",
    species_col="species"
)

# Per-isolate LSDD (requires df['species'])
df_lsdd = lsdd(
    df,
    lin_col="LINcode (C. jejuni / C. coli cgMLST v2)",
    species_col="species"
)
```

## Source-attribution colour scheme

Colours are hard-coded for consistency across figures:

- chicken: yellow
- ruminant: green
- pig: pink
- wild bird: purple
- other animal: grey
- human: near-black

See `linwalker/palette.py`.

## Notes on PubMLST exports

For PubMLST cgMLST exports, you usually do **not** need the full allele matrix for LINwalker.
A minimal export that includes `id`, `isolate`, `country`, `source`, `species`, and full `LINcode` is sufficient.

## License

MIT
