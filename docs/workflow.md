# LINwalker workflow (Campylobacter / PubMLST example)

This workflow uses *Campylobacter jejuni* / *C. coli* PubMLST exports (cgMLST v2 + LIN codes) as a worked example.
It is suitable for inclusion in the GitHub wiki.

---

## 0. Overview

LINwalker provides three related analyses:

1. **Diversification by scale**: unique LIN IDs vs LIN threshold (host-associated structure).
2. **Mixed-species LINs**: fraction of LIN clusters containing >1 species vs threshold (introgression signal).
3. **LSDD**: per-isolate **LIN Species Discordance Depth** (earliest LIN level where a lineage majority differs from isolate species).

---

## 1. Export from PubMLST

Recommended minimum fields (Provenance):
- id
- isolate
- country
- source
- species

Typing schemes:
- C. jejuni / C. coli cgMLST v2 (for LINcode)
- Select **all LIN code prefixes** (17/17) OR include full LINcode field
- (Optional) Export CAMP loci allele numbers if you want the cgMLST matrix for ML workflows

If you export CAMP loci, PubMLST output will be very wide; LINwalker can split it into smaller tables.

---

## 2. Clone + install

```bash
git clone https://github.com/Benizao1980/LINwalker
cd LINwalker

conda create -n linwalker python=3.11 -y
conda activate linwalker
pip install -r requirements.txt
```

---

## 3. Prepare the export (data cleaning + source binning)

```bash
python -m linwalker prep --input PATHSAFE_pubmlst_export.tsv.gz --outdir data/derived --prefix PATHSAFE

This collapses PubMLST `source` labels into stable bins: chicken, ruminant, pig, wild bird, human, other.
```

Outputs:
- `PATHSAFE_LINwalker_min.tsv` — minimal table for LINwalker analyses
- `PATHSAFE_metadata_only.tsv` — convenience metadata
- `PATHSAFE_cgMLST_matrix.tsv.gz` — CAMP loci matrix (only if present)

---

## 4. Diversification by host reservoir

```bash
python -m linwalker diversify   --input data/derived/PATHSAFE_LINwalker_min.tsv   --lin-col LINcode   --group-col source   --outdir results/diversification
```

Outputs:
- `diversification.tsv`
- `diversification.svg` / `diversification.png`

Interpretation:
- Divergence of curves across thresholds indicates scale-dependent host structure.
- Mid-range LIN levels often capture ecologically meaningful structure (not too coarse / not too fragmented).

---

## 5. Introgression summaries (mixed LINs + LSDD)

```bash
python -m linwalker introgress   --input data/derived/PATHSAFE_LINwalker_min.tsv   --lin-col LINcode   --species-col species   --outdir results/introgression
```

Outputs:
- `mixed_species.tsv` + plots
- `lsdd.tsv` + LSDD-by-source plots

Interpretation:
- Rising mixed-species fraction at intermediate LIN levels suggests structured introgression.
- Lower LSDD indicates deeper discordance (more hybrid-like at coarse scales).

---

## 6. Downstream integration

LINwalker outputs can be joined back to your metadata (or trees) by `id` / `isolate` and used to:
- justify attribution-relevant LIN thresholds
- interpret ML source predictions
- stratify introgression signals by reservoir, geography, or time
