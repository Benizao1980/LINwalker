# linwalker/stcc.py

from __future__ import annotations
import pandas as pd


def _lin_prefix(series: pd.Series, k: int) -> pd.Series:
    parts = series.astype(str).str.split("_")
    return parts.apply(lambda x: "_".join(x[:k]))


def stcc_concordance_by_threshold(
    df: pd.DataFrame,
    lin_col: str = "LINcode",
    st_col: str = "ST (MLST)",
    cc_col: str = "clonal_complex (MLST)",
    k_min: int = 1,
    k_max: int = 17,
    min_cluster_size: int = 2,
) -> pd.DataFrame:
    """
    Summarise concordance between LIN clusters and MLST ST / clonal complex across LIN thresholds.

    For each LIN level k:
      - form LIN clusters by k-prefix
      - compute per-cluster purity for ST and CC (max category proportion)
      - summarise (unweighted and weighted) plus proportions of "pure" clusters

    Parameters
    ----------
    min_cluster_size : int
        Exclude clusters smaller than this size from summaries (default 2).

    Returns
    -------
    DataFrame with one row per LIN_level and columns:
      - n_clusters
      - n_isolates_included
      - prop_pure_ST, prop_pure_CC
      - mean_purity_ST, mean_purity_CC
      - weighted_purity_ST, weighted_purity_CC
      - median_n_ST_per_cluster, median_n_CC_per_cluster
    """
    d = df.copy()
    for col in [lin_col, st_col, cc_col]:
        if col not in d.columns:
            raise ValueError(f"Missing required column: {col}")

    d[st_col] = d[st_col].astype(str).str.strip()
    d[cc_col] = d[cc_col].astype(str).str.strip()

    out_rows = []
    for k in range(k_min, k_max + 1):
        d["_lin_k"] = _lin_prefix(d[lin_col], k)
        g = d.groupby("_lin_k", dropna=False)

        # cluster sizes
        sizes = g.size().rename("cluster_size").reset_index()

        # per-cluster ST stats
        st_n = g[st_col].nunique().rename("n_ST")
        cc_n = g[cc_col].nunique().rename("n_CC")

        # purity = max frequency / size (ignore missing-like noticing: treat "nan" as category)
        st_purity = g[st_col].apply(lambda s: s.value_counts(dropna=False).iloc[0] / len(s)).rename("purity_ST")
        cc_purity = g[cc_col].apply(lambda s: s.value_counts(dropna=False).iloc[0] / len(s)).rename("purity_CC")

        per = pd.concat([st_n, cc_n, st_purity, cc_purity], axis=1).reset_index().merge(sizes, on="_lin_k")

        # filter small clusters
        per_f = per[per["cluster_size"] >= min_cluster_size].copy()
        if per_f.empty:
            out_rows.append({
                "LIN_level": k,
                "n_clusters": 0,
                "n_isolates_included": 0,
                "prop_pure_ST": 0.0,
                "prop_pure_CC": 0.0,
                "mean_purity_ST": 0.0,
                "mean_purity_CC": 0.0,
                "weighted_purity_ST": 0.0,
                "weighted_purity_CC": 0.0,
                "median_n_ST_per_cluster": 0.0,
                "median_n_CC_per_cluster": 0.0,
                "min_cluster_size": min_cluster_size,
            })
            continue

        n_clusters = int(per_f.shape[0])
        n_iso = int(per_f["cluster_size"].sum())

        prop_pure_ST = float((per_f["n_ST"] == 1).mean())
        prop_pure_CC = float((per_f["n_CC"] == 1).mean())

        mean_purity_ST = float(per_f["purity_ST"].mean())
        mean_purity_CC = float(per_f["purity_CC"].mean())

        weighted_purity_ST = float((per_f["purity_ST"] * per_f["cluster_size"]).sum() / n_iso)
        weighted_purity_CC = float((per_f["purity_CC"] * per_f["cluster_size"]).sum() / n_iso)

        out_rows.append({
            "LIN_level": k,
            "n_clusters": n_clusters,
            "n_isolates_included": n_iso,
            "prop_pure_ST": prop_pure_ST,
            "prop_pure_CC": prop_pure_CC,
            "mean_purity_ST": mean_purity_ST,
            "mean_purity_CC": mean_purity_CC,
            "weighted_purity_ST": weighted_purity_ST,
            "weighted_purity_CC": weighted_purity_CC,
            "median_n_ST_per_cluster": float(per_f["n_ST"].median()),
            "median_n_CC_per_cluster": float(per_f["n_CC"].median()),
            "min_cluster_size": min_cluster_size,
        })

    return pd.DataFrame(out_rows)
