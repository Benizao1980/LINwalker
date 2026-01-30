# linwalker/introgression.py

import pandas as pd
import numpy as np

def mixed_species_summary(df, lin_col, species_col):
    """
    Quantify mixed-species LIN clusters across thresholds.
    """
    df = df.copy()
    df["_lin_parts"] = df[lin_col].astype(str).str.split("_")
    max_k = df["_lin_parts"].apply(len).max()

    results = []

    for k in range(1, max_k + 1):
        df["_lin_k"] = df["_lin_parts"].apply(lambda x: tuple(x[:k]))
        grouped = df.groupby("_lin_k")[species_col]

        prop_mixed = (grouped.nunique() > 1).mean()

        entropy = grouped.apply(
            lambda s: -sum(
                p * np.log2(p)
                for p in s.value_counts(normalize=True)
            )
        ).mean()

        results.append({
            "LIN_level": k,
            "prop_mixed_species": prop_mixed,
            "mean_species_entropy": entropy
        })

    return pd.DataFrame(results)


def lsdd(df, lin_col, species_col):
    """
    Compute LIN Species Discordance Depth (LSDD) per isolate.
    """
    df = df.copy()
    df["_lin_parts"] = df[lin_col].astype(str).str.split("_")
    max_k = df["_lin_parts"].apply(len).max()

    lsdd_values = []

    for idx, row in df.iterrows():
        isolate_species = row[species_col]
        for k in range(1, max_k + 1):
            prefix = tuple(row["_lin_parts"][:k])
            cluster = df[df["_lin_parts"].apply(lambda x: tuple(x[:k]) == prefix)]
            majority_species = cluster[species_col].mode()[0]
            if majority_species != isolate_species:
                lsdd_values.append(k)
                break
        else:
            lsdd_values.append(max_k + 1)

    df["LSDD"] = lsdd_values
    return df
