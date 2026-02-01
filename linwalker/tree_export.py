from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import pandas as pd

from .utils import lin_prefix


@dataclass
class TreeExportResult:
    metadata: pd.DataFrame


def export_tree_metadata(
    df: pd.DataFrame,
    *,
    outdir: Path,
    lin_col: str = "lin_code",
    sample_col: str = "isolate",
    source_col: str = "source",
    species_col: str = "species",
    threshold: int = 12,
    extra_cols: Optional[str] = None,
) -> TreeExportResult:
    """Create a simple metadata table for iTOL/Microreact.

    The main field is a `LIN_<threshold>` column containing the LIN prefix,
    which can be used for colouring/labeling.
    """

    outdir.mkdir(parents=True, exist_ok=True)

    cols = [c for c in [sample_col, source_col, species_col, lin_col] if c in df.columns]
    m = df[cols].copy()
    m[f"LIN_{threshold}"] = m[lin_col].astype(str).apply(lambda s: lin_prefix(s, threshold))

    if extra_cols:
        for c in [x.strip() for x in extra_cols.split(",") if x.strip()]:
            if c in df.columns and c not in m.columns:
                m[c] = df[c]

    # Best-effort: set the sample id as the first column
    if sample_col in m.columns:
        m = m[[sample_col] + [c for c in m.columns if c != sample_col]]

    m.to_csv(outdir / "tree_metadata.tsv", sep="\t", index=False)
    return TreeExportResult(metadata=m)
