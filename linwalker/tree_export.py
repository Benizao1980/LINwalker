"""Exports for tree visualisation tools.

This module does not *draw* trees. Instead, it writes small metadata files that
can be imported into existing tree viewers.

Currently supported:
- Microreact: TSV with isolate id + columns to colour-by (e.g. source, species)
- iTOL: colourstrip file (simple, categorical)

You can colour an existing tree without re-rendering it in LINwalker.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import pandas as pd

from .palette import SOURCE_COLOURS


def export_microreact_tsv(
    df: pd.DataFrame,
    outpath: str | Path,
    id_col: str = "isolate",
    source_col: str = "source",
    species_col: str = "species",
    lin_col: str = "LINcode",
    extra_cols: Optional[list[str]] = None,
) -> Path:
    """Write a Microreact-compatible TSV.

    Microreact will import any TSV with a sample id column; users can then choose
    which column to colour by.
    """
    outpath = Path(outpath)
    outpath.parent.mkdir(parents=True, exist_ok=True)

    keep = [id_col, source_col, species_col, lin_col]
    if extra_cols:
        keep += [c for c in extra_cols if c in df.columns and c not in keep]

    out = df.loc[:, [c for c in keep if c in df.columns]].copy()
    out.to_csv(outpath, sep="\t", index=False)
    return outpath


def export_itol_colourstrip(
    df: pd.DataFrame,
    outpath: str | Path,
    id_col: str = "isolate",
    category_col: str = "source",
    title: str = "LINwalker",
    palette: Optional[dict[str, str]] = None,
) -> Path:
    """Write a minimal iTOL colourstrip file.

    iTOL colourstrip format:
      DATASET_COLORSTRIP
      SEPARATOR TAB
      DATASET_LABEL ...
      COLOR ...
      DATA
      <id> <hex> <label>

    This is intentionally simple and works well for source bins.
    """
    outpath = Path(outpath)
    outpath.parent.mkdir(parents=True, exist_ok=True)

    pal = palette or SOURCE_COLOURS

    tmp = df[[id_col, category_col]].copy()
    tmp[category_col] = tmp[category_col].astype(str).str.lower()
    tmp["hex"] = tmp[category_col].map(pal).fillna("#808080")

    lines = [
        "DATASET_COLORSTRIP",
        "SEPARATOR\tTAB",
        f"DATASET_LABEL\t{title}",
        "COLOR\t#000000",
        "DATA",
    ]
    for _, row in tmp.iterrows():
        lines.append(f"{row[id_col]}\t{row['hex']}\t{row[category_col]}")

    outpath.write_text("\n".join(lines) + "\n")
    return outpath
