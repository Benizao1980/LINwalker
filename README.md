# LINwalker

**LINwalker** is a lightweight Python toolkit for exploring hierarchical population structure and interspecies introgression using LIN codes (e.g., from PubMLST cgMLST schemes).

LINwalker treats LIN codes as ordered, multi-resolution descriptors rather than flat identifiers, enabling scale-aware analysis of:
- lineage diversification
- reservoir structure
- species boundary erosion (introgression)
- concordance with MLST ST and clonal complex

The package was developed and tested using Campylobacter jejuni / coli cgMLST data, but is applicable to any organism with LIN annotations.

## What it does

LINwalker provides five core analyses:

1. **Diversification by scale**
How the number of unique LIN clusters grows as LIN resolution increases (thresholds 1–17), stratified by source.

2. **Mixed-species LIN clusters**
Proportion of LIN clusters that contain more than one species at each LIN threshold — a scale-aware introgression signal.

3. **LSDD (LIN Species Discordance Depth)**
A per-isolate metric describing the earliest LIN level at which the isolate’s LIN cluster majority species differs from its assigned species.

4. **LIN <> MLST concordance**
Quantifies how well LIN clusters correspond to MLST sequence types (ST) and clonal complexes (CC) across LIN thresholds.

5. **Outbreak / public-health summaries**
   Descriptive outputs to support outbreak-style exploration: per-isolate LIN cluster
   size vs threshold (with optional boxplot+points), a table of the largest clusters
   at a chosen LIN threshold, and optional epi-curve by collection date.
   
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
conda install -c conda-forge \
  numpy \
  pandas \
  matplotlib \
  seaborn \
  scikit-learn
```

## Minimal plotting install (no Qt)

If you do not need interactive plotting:
```bash
conda install -c conda-forge matplotlib-base
```

This significantly reduces disk usage.

## Command-line usage

LINwalker ships with a small CLI. Run:

```bash
python -m linwalker --help
```

### 1) Prepare a PubMLST export (with source binning)

This step reconstructs a full `LINcode` column and collapses heterogeneous PubMLST `source` labels into stable bins:
**chicken, ruminant, pig, wild bird, human, other**.

```bash
python -m linwalker prep \
  --input PATHSAFE_pubmlst_export.tsv.gz \
  --outdir linwalker_run_v1_0_4/derived \
  --prefix PATHSAFE
```

This produces:
- `PATHSAFE_LINwalker_min.tsv` (LIN + metadata only)
- `PATHSAFE_metadata_only.tsv`
- `PATHSAFE_cgMLST_matrix.tsv.gz` (optional downstream use)

To preserve raw sources without binning:

```bash
python -m linwalker prep \
  --input PATHSAFE_pubmlst_export.tsv.gz \
  --outdir linwalker_run_v1_0_4/derived \
  --prefix PATHSAFE \
  --no-bin-sources
```

### 2) Diversification analysis

```bash
python -m linwalker diversify \
  --input linwalker_run_v1_0_4/derived/PATHSAFE_LINwalker_min.tsv \
  --outdir linwalker_run_v1_0_4/results_diversification
```

Outputs two versions by default:
- `diversification.png` / `.svg` Reservoir sources only (default ecological view)
- `diversification_all_sources.png` / `.svg` Includes human + other (complete epidemiological view)

### 3) Introgression analyses

```bash
python -m linwalker introgress \
  --input linwalker_run_v1_0_4/derived/PATHSAFE_LINwalker_min.tsv \
  --outdir linwalker_run_v1_0_4/results_introgression
```

Outputs:
- `mixed_species.png` / `.svg` Proportion of mixed-species LIN clusters vs LIN threshold (1–17)
- `lsdd_by_source.png / .svg` Distribution of LIN Species Discordance Depth by source

### 4) LIN <> MLST ST / clonal complex concordance

```bash
python -m linwalker stcc \
  --input linwalker_run_v1_0_4/derived/PATHSAFE_metadata_only.tsv \
  --outdir linwalker_run_v1_0_4/results_stcc

```

### 5) Outbreak / public-health summaries

Produces descriptive plots (cluster sizes vs LIN threshold, boxplots+points) and
tables of the largest LIN clusters at a chosen threshold. If you provide a date
column, LINwalker can also output a basic epi-curve.

```bash
python -m linwalker outbreak \
  --input linwalker_run_v1_0_4/derived/PATHSAFE_LINwalker_min.tsv \
  --outdir linwalker_run_v1_0_4/results_outbreak \
  --top-threshold 12 \
  --top-n 25 \
  --max-level 17 \
  --formats png svg
```

Outputs:
- `cluster_size_summary.*` median/IQR per-isolate cluster size vs LIN threshold
- `cluster_size_boxplot.*` boxplot + per-isolate jitter for cluster sizes
- `cluster_levels_summary.tsv` summary statistics per LIN level
- `top_clusters_t{threshold}.tsv` largest clusters at the chosen threshold
- `epi_curve_*.*` if you provide a `--date-col`

> Note: `stcc` expects columns named `ST (MLST)` and `clonal_complex (MLST)` if you feed it the `prep` output.

## Notes on interpretation (*Campylobacter*)
- LIN thresholds are strictly 1–17
- LIN codes are treated as cumulative prefixes, not independent columns
- Human isolates often dominate diversity curves and are therefore separated by default
- LSDD provides a scale-aware measure of species boundary erosion
- ST/CC concordance identifies LIN thresholds that approximate legacy typing units

### Interpreting LIN-based introgression metrics (LSDD)

LIN codes describe genetic relatedness at multiple nested scales, from very coarse (early LIN levels) to very fine (later LIN levels). This makes them useful not just for clustering isolates, but for asking where in the genetic hierarchy different biological signals appear.

#### *What is LSDD?*

LSDD (LIN Species Discordance Depth) is a per-isolate measure of how deep into the LIN hierarchy you have to go before species labels become inconsistent.

In practice, for each isolate:
1. LINwalker considers the isolate’s LIN clusters at each threshold (LIN 1 → LIN 17).
2. At each threshold, it asks:
    - *“What is the majority species among all isolates sharing this LIN prefix?”*
3. LSDD is defined as the earliest LIN level at which the isolate’s assigned species differs from the majority species of its LIN cluster.

If no discordance is observed at any level (LIN 1–17), the isolate is assigned the maximum value.

#### *How to interpret LSDD values*

LSDD is scale-aware by construction. Its biological interpretation depends on where discordance appears.

*Low LSDD values (early LIN levels)*
- Species discordance appears at coarse genetic scales
- Indicates deep or widespread mixing between species
- Suggests erosion of species boundaries that extends across broad lineages

In Campylobacter, this may reflect:
- long-term introgression
- shared ancestral structure
- extensive recombination across species boundaries

*High LSDD values (late LIN levels)*
- Species discordance only appears at fine genetic scales
- Most of the lineage structure remains species-consistent
- Mixing is localized or recent

This pattern is consistent with:
- occasional recombination events
- rare hybrid lineages
- spillover without sustained transmission

*Maximum LSDD (no discordance)*
- Species identity is consistent across all LIN levels
- No evidence of detectable interspecies mixing within the LIN hierarchy

#### *Why LSDD is useful*

Traditional approaches often treat introgression as a binary property (mixed vs not mixed). LSDD instead asks:
    *At what evolutionary scale does species mixing become visible?*

This allows you to distinguish between:
- deep, lineage-wide species boundary erosion
- versus shallow, fine-scale recombination events

Because LSDD is calculated per isolate, it can also be:
- summarised by source or host
- compared across datasets
- linked to other metadata (e.g. ecology, geography, clinical status)

#### *Important caveats*
- LSDD does not identify specific recombination tracts or donor lineages.
- It should be interpreted as a population-structural signal, not direct mechanistic evidence.
- Values depend on the resolution and composition of the reference dataset.

LSDD is therefore best used as:
- a comparative, scale-aware summary of species boundary stability within a dataset.

## Source-attribution colour scheme

Colours are hard-coded for consistency across figures:

- chicken: yellow
- ruminant: green
- pig: pink
- wild bird: purple
- other: grey
- human: near-black

See `linwalker/palette.py`.

## Citation
Please cite Parfitt et al. (*In preparation*).
