# linwalker/diversification.py

from __future__ import annotations
import pandas as pd


def lin_diversification(df: pd.DataFrame, lin_col: str, group_col: str) -> pd.DataFrame:
    """
    Count unique LIN identifiers across hierarchical thresholds.

    Parameters
    ----------
    df : DataFrame
        Input dataframe containing LIN codes and grouping labels.
    lin_col : str
        Column containing full LIN code as underscore-separated levels.
    group_col : str
        Column defining groups (e.g. source/reservoir).

    Returns
    -------
    DataFrame with columns:
        LIN_level, group, n_unique_LINs
    """
    d = df.copy()
    d["_lin_parts"] = d[lin_col].astype(str).str.split("_")
    max_k = int(d["_lin_parts"].apply(len).max())

    rows = []
    for k in range(1, max_k + 1):
        d["_lin_k"] = d["_lin_parts"].apply(lambda x: tuple(x[:k]))
        for grp, sub in d.groupby(group_col):
            rows.append({
                "LIN_level": k,
                "group": grp,
                "n_unique_LINs": int(sub["_lin_k"].nunique())
            })

    out = pd.DataFrame(rows).sort_values(["group", "LIN_level"]).reset_index(drop=True)
    return out
