"""LIN <> MLST concordance.

At each LIN threshold *k*, we form LIN clusters using the LIN prefix of length *k*.
We then compare those clusters to MLST sequence type (ST) and clonal complex (CC).

Outputs (per threshold):
- purity_ST / purity_CC: proportion of LIN clusters that are *pure* (single ST/CC)
- weighted_purity_ST / weighted_purity_CC: purity weighted by cluster sizes
- ARI_ST / ARI_CC: Adjusted Rand Index comparing LIN clusters to ST/CC labels

Notes
-----
* Purity answers: "are LIN clusters internally consistent?"
* ARI answers: "do LIN clusters reproduce the same partition as ST/CC?" (up to relabelling)

This is a descriptive mapping, not a claim that one scheme is "correct".
"""

from __future__ import annotations

from typing import Iterable, List, Optional, Tuple

import numpy as np
import pandas as pd
from sklearn.metrics import adjusted_rand_score

from .utils import infer_max_level_from_lincode, parse_thresholds


def lin_prefix(code: str, k: int) -> str:
    parts = str(code).split("_")
    return "_".join(parts[:k])


def _resolve_column(df: pd.DataFrame, preferred: str, alternatives: List[str]) -> str:
    """Return the first matching column name.

    We accept common PubMLST naming variants because exports differ.
    """
    if preferred in df.columns:
        return preferred
    for alt in alternatives:
        if alt in df.columns:
            return alt
    # also try a case-insensitive match
    lower = {c.lower(): c for c in df.columns}
    if preferred.lower() in lower:
        return lower[preferred.lower()]
    for alt in alternatives:
        if alt.lower() in lower:
            return lower[alt.lower()]
    raise ValueError(
        f"Missing required column: {preferred}. Available columns: {', '.join(df.columns)}"
    )


def _cluster_purity(labels: pd.Series) -> float:
    # labels: one column within a cluster
    vals = labels.dropna().astype(str)
    if len(vals) == 0:
        return 0.0
    return 1.0 if vals.nunique() == 1 else 0.0


def _weighted_purity(cluster_sizes: np.ndarray, purities: np.ndarray) -> float:
    if cluster_sizes.sum() == 0:
        return 0.0
    return float((cluster_sizes * purities).sum() / cluster_sizes.sum())


def stcc_concordance_by_threshold(
    df: pd.DataFrame,
    lin_col: str = "LINcode",
    st_col: str = "ST (MLST)",
    cc_col: str = "clonal_complex (MLST)",
    thresholds: Optional[Iterable[int]] = None,
    min_cluster_size: int = 1,
) -> pd.DataFrame:
    """Compute LIN<>ST/CC concordance summaries across thresholds."""

    if lin_col not in df.columns:
        raise ValueError(f"Missing required LIN column: {lin_col}")

    # Be permissive about ST/CC column names (common PubMLST exports)
    st_col = _resolve_column(df, st_col, ["ST", "mlst_st", "MLST_ST", "ST (MLST)"])
    cc_col = _resolve_column(
        df,
        cc_col,
        ["clonal_complex", "CC", "mlst_cc", "clonal_complex (MLST)", "CC (MLST)"],
    )

    max_level = infer_max_level_from_lincode(df[lin_col])
    ks = parse_thresholds(thresholds, max_level=max_level)

    # Pre-coerce to strings (stable ARI behaviour)
    st_labels = df[st_col].astype(str).fillna("NA")
    cc_labels = df[cc_col].astype(str).fillna("NA")

    rows = []
    for k in ks:
        clusters = df[lin_col].astype(str).apply(lambda s: lin_prefix(s, k))

        tmp = df[[lin_col]].copy()
        tmp["cluster"] = clusters
        tmp["ST"] = st_labels
        tmp["CC"] = cc_labels

        # Filter very small clusters if requested
        if min_cluster_size > 1:
            sizes = tmp["cluster"].value_counts()
            keep = sizes[sizes >= min_cluster_size].index
            tmp = tmp[tmp["cluster"].isin(keep)].copy()

        if tmp.empty:
            rows.append(
                {
                    "threshold": k,
                    "n_clusters": 0,
                    "purity_ST": 0.0,
                    "purity_CC": 0.0,
                    "weighted_purity_ST": 0.0,
                    "weighted_purity_CC": 0.0,
                    "ARI_ST": np.nan,
                    "ARI_CC": np.nan,
                }
            )
            continue

        grouped = tmp.groupby("cluster", sort=False)
        cluster_sizes = grouped.size().to_numpy()

        purity_st = grouped["ST"].apply(_cluster_purity).to_numpy()
        purity_cc = grouped["CC"].apply(_cluster_purity).to_numpy()

        # ARI compares partitions of isolates, not cluster-internal stats.
        ari_st = adjusted_rand_score(tmp["ST"], tmp["cluster"])
        ari_cc = adjusted_rand_score(tmp["CC"], tmp["cluster"])

        rows.append(
            {
                "threshold": k,
                "n_clusters": int(grouped.ngroups),
                "purity_ST": float(purity_st.mean()),
                "purity_CC": float(purity_cc.mean()),
                "weighted_purity_ST": _weighted_purity(cluster_sizes, purity_st),
                "weighted_purity_CC": _weighted_purity(cluster_sizes, purity_cc),
                "ARI_ST": float(ari_st),
                "ARI_CC": float(ari_cc),
            }
        )

    return pd.DataFrame(rows)
