from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

from .utils import coerce_source, lin_prefix, parse_thresholds


@dataclass
class DiversificationResult:
    table: pd.DataFrame


@dataclass
class RarefactionResult:
    """Rarefied diversification summary.

    table columns:
      - source
      - lin_level
      - n_unique_mean
      - n_unique_sd
      - n_per_group
      - n_reps
    """

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


def _rarefy_one(
    df: pd.DataFrame,
    group_col: str,
    n_per_group: int,
    seed: int,
) -> pd.DataFrame:
    """Return a rarefied dataframe with n_per_group rows per group."""
    rng = np.random.default_rng(seed)
    parts = []
    for g, x in df.groupby(group_col, dropna=False):
        if len(x) < n_per_group:
            continue
        idx = rng.choice(x.index.to_numpy(), size=n_per_group, replace=False)
        parts.append(x.loc[idx])
    if not parts:
        return df.iloc[0:0].copy()
    return pd.concat(parts, axis=0)


def lin_diversification_rarefied(
    df: pd.DataFrame,
    lin_col: str = "LINcode",
    group_col: str = "source",
    thresholds: str = "1-17",
    max_level: int = 17,
    n_per_group: Optional[int] = None,
    n_reps: int = 100,
    seed: int = 123,
) -> RarefactionResult:
    """Rarefied diversification curves.

    For each group, subsample an equal number of isolates (n_per_group) and
    compute the number of unique LIN prefixes across thresholds. Repeat
    n_reps times and report mean±sd.

    If n_per_group is not provided, uses the minimum group size across groups.
    """

    if lin_col not in df.columns:
        raise ValueError(f"Missing LIN column: {lin_col}")
    if group_col not in df.columns:
        raise ValueError(f"Missing group column: {group_col}")

    from .utils import lin_prefix, parse_thresholds

    work = df[[lin_col, group_col]].copy()
    work = work.dropna(subset=[lin_col, group_col])

    # Determine rarefaction depth
    sizes = work.groupby(group_col).size()
    if sizes.empty:
        raise ValueError("No rows available after filtering.")
    if n_per_group is None:
        n_per_group = int(sizes.min())
    if n_per_group <= 0:
        raise ValueError("n_per_group must be >0")

    ks = parse_thresholds(thresholds, max_level)

    # Collect replicate results
    rows: List[Dict] = []
    for r in range(n_reps):
        sub = _rarefy_one(work, group_col=group_col, n_per_group=n_per_group, seed=seed + r)
        if sub.empty:
            continue
        for k in ks:
            tmp = sub.copy()
            tmp["lin_prefix"] = tmp[lin_col].astype(str).map(lambda s: lin_prefix(s, k))
            grp = tmp.groupby(group_col)["lin_prefix"].nunique(dropna=True)
            for g, nuniq in grp.items():
                rows.append({"rep": r, "source": g, "lin_level": k, "n_unique": int(nuniq)})

    res = pd.DataFrame(rows)
    if res.empty:
        raise ValueError("Rarefaction produced no results (check group sizes / n_per_group).")

    summ = (
        res.groupby(["source", "lin_level"], as_index=False)["n_unique"]
        .agg([("n_unique_mean", "mean"), ("n_unique_sd", "std")])
        .reset_index()
    )
    summ["n_per_group"] = n_per_group
    summ["n_reps"] = n_reps

    # Ensure stable ordering
    summ = summ.sort_values(["source", "lin_level"]).reset_index(drop=True)
    return RarefactionResult(table=summ)
