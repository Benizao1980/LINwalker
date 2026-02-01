from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional

import pandas as pd

from .utils import coerce_source, lin_prefix, parse_thresholds


@dataclass
class DiversificationResult:
    table: pd.DataFrame


def lin_diversification(
    df: pd.DataFrame,
    *,
    lin_col: str = "lin_code",
    group_col: str = "source",
    thresholds: Optional[str] = None,
    max_level: int = 17,
) -> DiversificationResult:
    """Compute number of unique LIN prefixes per LIN level, optionally per group."""

    if lin_col not in df.columns:
        raise ValueError(f"Missing LIN column: {lin_col}")
    if group_col not in df.columns:
        raise ValueError(f"Missing group column: {group_col}")

    ks = parse_thresholds(thresholds, max_level)

    # normalize
    work = df[[lin_col, group_col]].copy()
    work[group_col] = work[group_col].apply(coerce_source)

    records: List[Dict] = []
    for k in ks:
        tmp = work.copy()
        tmp["lin_prefix"] = tmp[lin_col].astype(str).map(lambda s: lin_prefix(s, k))
        # unique per group
        grp = tmp.groupby(group_col)["lin_prefix"].nunique(dropna=True).reset_index()
        grp.insert(0, "lin_level", k)
        grp.rename(columns={"lin_prefix": "unique_lin"}, inplace=True)
        records.append(grp)

    out = pd.concat(records, ignore_index=True)

    # Standardise column names for downstream plotting/CLI.
    if group_col != "source":
        out = out.rename(columns={group_col: "source"})
    if "unique_lin" in out.columns:
        out = out.rename(columns={"unique_lin": "n_unique"})
    out = out.sort_values(["source", "lin_level"])

    return DiversificationResult(table=out)
