from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from .utils import coerce_source


@dataclass
class PrepResult:
    minimal: pd.DataFrame
    metadata_only: pd.DataFrame
    cgmlst_matrix: pd.DataFrame | None = None


def prep_pubmlst_export(
    input_path: Path,
    *,
    outdir: Path,
    prefix: str = "RUN",
    bin_sources: bool = True,
    lin_col: str = "LINcode",
    cgst_col: str = "cgST",
    st_col: str = "ST",
    cc_col: str = "clonal_complex",
    species_col: str = "species",
    source_col: str = "source",
    country_col: str = "country",
    date_col: str = "date",
    sample_col: str = "isolate",
    keep_cgmlst_matrix: bool = True,
) -> PrepResult:
    """Prepare a PubMLST export for LINwalker.

    Outputs are written directly under *outdir*. The function standardises common
    metadata column names, optionally bins source labels, and can retain the wide
    CAMP-locus matrix separately.
    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    input_path = Path(input_path)
    # PubMLST exports used by LINwalker are normally TSV/TSV.GZ.
    sep = "," if input_path.name.lower().endswith(".csv") else "	"
    df = pd.read_csv(input_path, sep=sep, low_memory=False)

    rename = {}
    for src, dst in [
        (sample_col, "id"),
        (lin_col, "lin_code"),
        (cgst_col, "cgST"),
        (st_col, "ST"),
        (cc_col, "clonal_complex"),
        (species_col, "species"),
        (source_col, "source"),
        (country_col, "country"),
        (date_col, "collection_date"),
    ]:
        if src in df.columns and src != dst:
            rename[src] = dst
    if rename:
        df = df.rename(columns=rename)

    if "id" not in df.columns:
        for alt in ["isolate", "isolate_id", "isolateID", "name"]:
            if alt in df.columns:
                df = df.rename(columns={alt: "id"})
                break
    if "id" not in df.columns:
        df.insert(0, "id", [f"sample_{i+1}" for i in range(len(df))])

    if "source" in df.columns and bin_sources:
        df["source"] = df["source"].map(coerce_source)

    if "collection_date" in df.columns:
        df["collection_date"] = pd.to_datetime(df["collection_date"], errors="coerce")

    minimal_cols = [
        c for c in [
            "id", "lin_code", "species", "source", "country",
            "collection_date", "ST", "clonal_complex", "cgST",
        ]
        if c in df.columns
    ]
    min_df = df[minimal_cols].copy()

    meta_cols = [
        c for c in [
            "id", "species", "source", "country", "collection_date",
            "ST", "clonal_complex", "cgST",
        ]
        if c in min_df.columns
    ]
    meta_df = min_df[meta_cols].copy()

    min_out = outdir / f"{prefix}_LINwalker_min.tsv"
    meta_out = outdir / f"{prefix}_metadata_only.tsv"
    min_df.to_csv(min_out, sep="	", index=False)
    meta_df.to_csv(meta_out, sep="	", index=False)

    cg_df = None
    locus_cols = [c for c in df.columns if str(c).startswith("CAMP")]
    if keep_cgmlst_matrix and locus_cols:
        cg_df = df[["id"] + locus_cols].copy()
        cg_df.to_csv(
            outdir / f"{prefix}_cgMLST_matrix.tsv.gz",
            sep="	",
            index=False,
            compression="gzip",
        )

    return PrepResult(
        minimal=min_df,
        metadata_only=meta_df,
        cgmlst_matrix=cg_df,
    )
