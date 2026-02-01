from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional

import pandas as pd

from .utils import coerce_source


@dataclass
class PrepResult:
    minimal: pd.DataFrame
    metadata_only: pd.DataFrame


def prep_pubmlst_export(
    input_path: Path,
    *,
    outdir: Path,
    lin_col: str = "LINcode",
    cgst_col: str = "cgST",
    st_col: str = "ST",
    cc_col: str = "clonal_complex",
    species_col: str = "species",
    source_col: str = "source",
    country_col: str = "country",
    date_col: str = "collection_date",
) -> PrepResult:
    """Prepare a PubMLST export table for downstream LINwalker analyses.

    Writes:
      - derived/PATHSAFE_LINwalker_min.tsv
      - derived/PATHSAFE_metadata_only.tsv

    Returns DataFrames in memory too.
    """

    outdir.mkdir(parents=True, exist_ok=True)
    derived = outdir / "derived"
    derived.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(input_path, sep="\t", low_memory=False)

    rename = {
        lin_col: "lin_code",
        cgst_col: "cgST",
        st_col: "ST",
        cc_col: "clonal_complex",
        species_col: "species",
        source_col: "source",
        country_col: "country",
        date_col: "collection_date",
    }
    for k, v in rename.items():
        if k in df.columns:
            df = df.rename(columns={k: v})

    # Normalise key columns
    if "source" in df.columns:
        df["source"] = df["source"].map(coerce_source)

    # Parse dates if present
    if "collection_date" in df.columns:
        df["collection_date"] = pd.to_datetime(df["collection_date"], errors="coerce")

    keep_min = [c for c in ["id", "lin_code", "species", "source", "country", "collection_date", "ST", "clonal_complex", "cgST"] if c in df.columns]
    if "id" not in df.columns:
        # Try a couple of common PubMLST columns
        for alt in ["isolate", "isolate_id", "isolateID", "name"]:
            if alt in df.columns:
                df = df.rename(columns={alt: "id"})
                break
        if "id" not in df.columns:
            df.insert(0, "id", [f"sample_{i+1}" for i in range(len(df))])
            keep_min = ["id"] + keep_min

    min_df = df.loc[:, dict.fromkeys(keep_min)].copy()
    meta_cols = [c for c in ["id", "species", "source", "country", "collection_date", "ST", "clonal_complex", "cgST"] if c in min_df.columns]
    meta_df = min_df.loc[:, meta_cols].copy()

    min_out = derived / "PATHSAFE_LINwalker_min.tsv"
    meta_out = derived / "PATHSAFE_metadata_only.tsv"
    min_df.to_csv(min_out, sep="\t", index=False)
    meta_df.to_csv(meta_out, sep="\t", index=False)

    return PrepResult(minimal=min_df, metadata_only=meta_df)
