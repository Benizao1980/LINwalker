# linwalker/introgression.py

from __future__ import annotations
import pandas as pd


def mixed_species_summary(df: pd.DataFrame, lin_col: str, species_col: str) -> pd.DataFrame:
    """
    Quantify mixed-species LIN clusters across thresholds.

    Output columns:
      - LIN_level
      - prop_mixed_species
      - n_clusters
      - n_mixed_clusters
    """
    d = df.copy()
    d["_lin_parts"] = d[lin_col].astype(str).str.split("_")
    max_k = int(d["_lin_parts"].apply(len).max())

    out = []
    for k in range(1, max_k + 1):
        keys = d["_lin_parts"].apply(lambda x: "_".join(x[:k]))
        tmp = pd.DataFrame({"key": keys, "species": d[species_col].astype(str)})
        nunq = tmp.groupby("key")["species"].nunique()
        n_clusters = int(nunq.shape[0])
        n_mixed = int((nunq > 1).sum())
        out.append({
            "LIN_level": k,
            "prop_mixed_species": (n_mixed / n_clusters) if n_clusters else 0.0,
            "n_clusters": n_clusters,
            "n_mixed_clusters": n_mixed
        })

    return pd.DataFrame(out)


def lsdd(df: pd.DataFrame, lin_col: str, species_col: str) -> pd.DataFrame:
    """
    Compute LIN Species Discordance Depth (LSDD) per isolate.

    LSDD is the earliest LIN level where the isolate's LIN cluster majority species
    differs from the isolate's assigned species. If no discordance occurs across
    all levels, LSDD = max_level + 1.
    """
    d = df.copy()
    d["_lin_parts"] = d[lin_col].astype(str).str.split("_")
    max_k = int(d["_lin_parts"].apply(len).max())

    parts = d["_lin_parts"].tolist()
    species = d[species_col].astype(str).tolist()

    # Precompute prefix keys for each k
    prefix_keys = []
    for k in range(1, max_k + 1):
        prefix_keys.append(["_".join(p[:k]) for p in parts])

    # Majority species per prefix at each k
    majority_by_k = []
    for k in range(1, max_k + 1):
        tmp = pd.DataFrame({"key": prefix_keys[k-1], "species": species})
        maj = tmp.groupby("key")["species"].agg(lambda s: s.mode().iloc[0])
        majority_by_k.append(maj)

    lsdd_vals = []
    for i in range(len(parts)):
        sp_i = species[i]
        found = False
        for k in range(1, max_k + 1):
            key = prefix_keys[k-1][i]
            if majority_by_k[k-1].loc[key] != sp_i:
                lsdd_vals.append(k)
                found = True
                break
        if not found:
            lsdd_vals.append(max_k + 1)

    d["LSDD"] = lsdd_vals
    return d
