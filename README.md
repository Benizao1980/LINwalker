# LINwalker

**LINwalker** is a lightweight Python toolkit for exploring hierarchical population structure and interspecies introgression using LIN codes (e.g., from PubMLST cgMLST schemes).

LINwalker treats LIN codes as **ordered, multi-resolution descriptors** rather than flat identifiers, enabling scale-aware analysis of host association, lineage structure, and recombination/introgression.

Repository: https://github.com/Benizao1980/LINwalker

## What it does

LINwalker provides three core analyses:

1. **Diversification by scale** — how the number of unique LIN IDs grows as you increase LIN resolution.
2. **Mixed-species LINs** — how often LIN prefixes contain more than one species (scale-aware introgression signal).
3. **LSDD (LIN Species Discordance Depth)** — per-isolate metric for the earliest LIN level at which its LIN cluster’s majority species differs from its assigned species.

## Getting started (clone → install → run)

```bash
git clone https://github.com/Benizao1980/LINwalker
cd LINwalker
```

### Install dependencies

**Conda (recommended):**
```bash
conda create -n linwalker python=3.11 -y
conda activate linwalker
pip install -r requirements.txt
```

**Or pip/venv:**
```bash
python -m venv .venv
source .venv/bin/activate   # macOS/Linux
# .venv\Scripts\activate    # Windows PowerShell
pip install -r requirements.txt
```

## Command-line usage

LINwalker ships with a small CLI. Run:

```bash
python -m linwalker --help
```

### 1) Prepare a PubMLST export (with source binning)

This step reconstructs a full `LINcode` column and collapses heterogeneous PubMLST `source` labels into stable bins:
**chicken, ruminant, pig, wild bird, human, other**.

```bash
python -m linwalker prep --input PATHSAFE_pubmlst_export.tsv.gz --outdir data/derived --prefix PATHSAFE
```

To preserve raw sources without binning:

```bash
python -m linwalker prep --input PATHSAFE_pubmlst_export.tsv.gz --outdir data/derived --prefix PATHSAFE --no-bin-sources
```

### 2) Diversification plot

```bash
python -m linwalker diversify --input data/derived/PATHSAFE_LINwalker_min.tsv --lin-col LINcode --group-col source --outdir results/diversification
```

### 3) Introgression summaries

```bash
python -m linwalker introgress --input data/derived/PATHSAFE_LINwalker_min.tsv --lin-col LINcode --species-col species --outdir results/introgression
```

## Python usage (API)

```python
from linwalker.prep import prepare_pubmlst_export
from linwalker.diversification import lin_diversification
from linwalker.introgression import mixed_species_summary, lsdd

lin_df, meta_df, cg_df = prepare_pubmlst_export("PATHSAFE_pubmlst_export.tsv.gz")

div = lin_diversification(lin_df, lin_col="LINcode", group_col="source")
mix = mixed_species_summary(lin_df, lin_col="LINcode", species_col="species")
lin_df_lsdd = lsdd(lin_df, lin_col="LINcode", species_col="species")
```

## Source-attribution colour scheme

Colours are hard-coded for consistency across figures:

- chicken: yellow
- ruminant: green
- pig: pink
- wild bird: purple
- other: grey
- human: near-black

See `linwalker/palette.py`.

## License

MIT
