# linwalker/diversification.py

import pandas as pd

def lin_diversification(df, lin_col, group_col):
    """
    Count unique LIN identifiers across hierarchical thresholds.

    Returns a dataframe with columns:
    - LIN_level
    - group
    - n_unique_LINs
    """
    df = df.copy()
    df["_lin_parts"] = df[lin_col].astype(str).str.split("_")
    max_k = df["_lin_parts"].apply(len).max()

    records = []

    for k in range(1, max_k + 1):
        df["_lin_k"] = df["_lin_parts"].apply(lambda x: tuple(x[:k]))
        for grp, sub in df.groupby(group_col):
            records.append({
                "LIN_level": k,
                "group": grp,
                "n_unique_LINs": sub["_lin_k"].nunique()
            })

    return pd.DataFrame(records)
