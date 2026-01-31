# linwalker/prep.py
"""
Data preparation utilities for LINwalker.

Primary use-case: cleaning PubMLST exports (TSV/TSV.GZ) into small, analysis-ready
tables suitable for LIN-based structure and introgression analyses.

Key features
- Reconstructs full LINcode from PubMLST LINcode prefix columns (1..17) if needed
- Handles PubMLST exports where LINcode[n] columns are *cumulative* (each cell already
  contains an underscore-separated code up to depth n)
- Normalises species labels
- Bins noisy source labels into stable source categories for plotting and attribution:
  chicken, ruminant, pig, wild bird, human, other
"""

from __future__ import annotations

import gzip
import re
from pathlib import Path
from typing import Optional, Tuple

import pandas as pd


def _open_text(path: str):
    """Open plain text or gzipped text file."""
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", errors="replace")
    return open(path, "rt", errors="replace")


def _normalise_text_series(s: pd.Series) -> pd.Series:
    s = s.astype(str).str.lower().str.strip()
    s = s.replace({"nan": pd.NA, "none": pd.NA, "null": pd.NA, "": pd.NA})
    return s


def bin_source_labels(source: pd.Series) -> pd.Series:
    """
    Collapse heterogeneous PubMLST 'source' labels into stable categories.

    Output categories:
      - chicken
      - ruminant
      - pig
      - wild bird
      - human
      - other
    """
    s = _normalise_text_series(source).fillna("")
    out = pd.Series(["other"] * len(s), index=s.index, dtype="object")

    # Human
    human_pat = r"(?:human|patient|stool|faec|fec|diarr|clinic|hospital|case)"
    out[s.str.contains(human_pat, regex=True)] = "human"

    # Chicken / poultry
    chicken_pat = r"(?:chicken|poultry|broiler|layer|hen|rooster|turkey)"
    out[s.str.contains(chicken_pat, regex=True)] = "chicken"

    # Pig
    pig_pat = r"(?:pig|swine|porcine)"
    out[s.str.contains(pig_pat, regex=True)] = "pig"

    # Ruminant (cattle, sheep, goats, etc.)
    ruminant_pat = r"(?:ruminant|cattle|cow|bovine|beef|calf|sheep|goat|caprine|ovine|lamb|dairy)"
    out[s.str.contains(ruminant_pat, regex=True)] = "ruminant"

    # Wild birds (explicitly wild)
    wildbird_pat = r"(?:wild\s*bird|wildbird|gull|crow|sparrow|pigeon|duck|goose|seabird|wader)"
    out[s.str.contains(wildbird_pat, regex=True)] = "wild bird"

    return out


def prepare_pubmlst_export(
    path: str,
    outdir: Optional[str] = None,
    prefix: Optional[str] = None,
    source_col: str = "source",
    species_col: str = "species",
    bin_sources: bool = True,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Prepare a PubMLST export for LINwalker analyses.

    Returns (lin_df, meta_df, cg_df).

    Notes
    -----
    PubMLST may export LIN prefixes as either:
      - "LINcode[1]" ... "LINcode[17]"
      - "LINcode[1] (scheme)" ... etc

    Critically, PubMLST prefix columns are often *cumulative* strings:
      LINcode[3] cell == "0_1_9"
    In that case, the deepest column (e.g. LINcode[17]) is already the full LINcode.
    """
    path = str(path)
    with _open_text(path) as f:
        df = pd.read_csv(f, sep="\t", low_memory=False)

    # Normalise species/source (if present)
    if species_col in df.columns:
        df[species_col] = _normalise_text_series(df[species_col])
    if source_col in df.columns:
        df[source_col] = _normalise_text_series(df[source_col])

    # Bin sources into stable categories (optional)
    if bin_sources and source_col in df.columns:
        df["source_raw"] = df[source_col]
        df[source_col] = bin_source_labels(df[source_col])

    cols = list(df.columns)

    # Identify LIN columns
    full_lin = [c for c in cols if c.startswith("LINcode (")]
    prefix_lin = [c for c in cols if re.match(r"^LINcode\[\d+\]", c)]  # allows suffix "(scheme)"

    if full_lin:
        lin_col = full_lin[0]
        df["LINcode"] = df[lin_col].astype(str)
    elif prefix_lin:
        prefix_lin = sorted(prefix_lin, key=lambda x: int(re.findall(r"\d+", x)[0]))
        deepest = prefix_lin[-1]
        sample = str(df[deepest].iloc[0])
        if "_" in sample:
            # cumulative prefixes -> deepest already full code
            df["LINcode"] = df[deepest].astype(str)
        else:
            # single-level tokens -> join
            df["LINcode"] = df[prefix_lin].astype(str).agg("_".join, axis=1)
    elif "LINcode" in cols:
        df["LINcode"] = df["LINcode"].astype(str)
    else:
        raise ValueError("No LINcode column found. Expected 'LINcode (...)' or 'LINcode[n]' prefix columns.")

    # Minimal linwalker table
    keep_min = [c for c in ["id", "isolate", "country", source_col, species_col, "LINcode"] if c in df.columns]
    lin_df = df[keep_min].copy()

    # Metadata-only table
    keep_meta = [
        "id", "isolate", "country", source_col, species_col,
        "ST (MLST)", "clonal_complex (MLST)",
        "cgST (C. jejuni / C. coli cgMLST v2)",
        "LINcode"
    ]
    if "source_raw" in df.columns:
        keep_meta.insert(4, "source_raw")
    meta_df = df[[c for c in keep_meta if c in df.columns]].copy()

    # cgMLST matrix (CAMP loci)
    camp_cols = [c for c in cols if re.match(r"^CAMP\d{4}$", c)]
    if camp_cols:
        cg_df = df[[c for c in ["id", "isolate", source_col] if c in df.columns] + camp_cols].copy()
    else:
        cg_df = pd.DataFrame()

    # Write outputs if requested
    if outdir:
        outdir_p = Path(outdir)
        outdir_p.mkdir(parents=True, exist_ok=True)
        if prefix is None:
            base = Path(path).name
            prefix = re.sub(r"\.(tsv|txt)(\.gz)?$", "", base, flags=re.IGNORECASE)

        lin_path = outdir_p / f"{prefix}_LINwalker_min.tsv"
        meta_path = outdir_p / f"{prefix}_metadata_only.tsv"
        lin_df.to_csv(lin_path, sep="\t", index=False)
        meta_df.to_csv(meta_path, sep="\t", index=False)

        if not cg_df.empty:
            cg_path = outdir_p / f"{prefix}_cgMLST_matrix.tsv.gz"
            with gzip.open(cg_path, "wt") as gz:
                cg_df.to_csv(gz, sep="\t", index=False)

    return lin_df, meta_df, cg_df
