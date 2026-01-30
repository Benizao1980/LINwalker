# linwalker/introgression.py

import pandas as pd
import numpy as np

def mixed_species_summary(df: pd.DataFrame, lin_col: str, species_col: str) -> pd.DataFrame:
    """Quantify mixed-species LIN clusters across thresholds.

    Returns a dataframe with:
        LIN_level, prop_mixed_species, mean_species_entropy
    """
    d = df.copy()
    d["_lin_parts"] = d[lin_col].astype(str).str.split("_")
    max_k = d["_lin_parts"].apply(len).max()

    out = []
    for k in range(1, max_k + 1):
        d["_lin_k"] = d["_lin_parts"].apply(lambda x: tuple(x[:k]))
        grouped = d.groupby("_lin_k")[species_col]

        prop_mixed = (grouped.nunique() > 1).mean()

        def entropy(s):
            p = s.value_counts(normalize=True)
            return float(-(p * np.log2(p)).sum()) if len(p) else 0.0

        mean_entropy = grouped.apply(entropy).mean()

        out.append({
            "LIN_level": k,
            "prop_mixed_species": float(prop_mixed),
            "mean_species_entropy": float(mean_entropy)
        })

    return pd.DataFrame(out)


def lsdd(df: pd.DataFrame, lin_col: str, species_col: str) -> pd.DataFrame:
    """Compute LIN Species Discordance Depth (LSDD) per isolate.

    LSDD is the earliest LIN level where the isolate's LIN cluster majority species
    differs from the isolate's assigned species. If no discordance occurs across
    all levels, LSDD = max_level + 1.
    """
    d = df.copy()
    d["_lin_parts"] = d[lin_col].astype(str).str.split("_")
    max_k = d["_lin_parts"].apply(len).max()

    # Precompute cluster keys for speed
    lin_parts = d["_lin_parts"].tolist()
    species = d[species_col].tolist()

    lsdd_vals = []
    for i, (parts_i, sp_i) in enumerate(zip(lin_parts, species)):
        found = False
        for k in range(1, max_k + 1):
            prefix = tuple(parts_i[:k])
            # cluster membership
            idxs = [j for j, parts_j in enumerate(lin_parts) if tuple(parts_j[:k]) == prefix]
            cluster_species = [species[j] for j in idxs]
            # majority species (mode)
            maj = pd.Series(cluster_species).mode().iloc[0]
            if maj != sp_i:
                lsdd_vals.append(k)
                found = True
                break
        if not found:
            lsdd_vals.append(max_k + 1)

    d["LSDD"] = lsdd_vals
    return d
