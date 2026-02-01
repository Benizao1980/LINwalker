"""Diversification by LIN threshold.

For each LIN level k, we define a LIN cluster as the prefix of length k.
We then count unique clusters either overall or within each source group.

The output is a long table with columns:
- threshold (int)
- source (str)
- n_clusters (int)

Plots are saved as PNG + SVG by default.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import pandas as pd

from .plotting import plot_diversification_curve
from .utils import infer_max_level_from_lincode, parse_thresholds


def lin_diversification(
    df: pd.DataFrame,
    lin_col: str = "LINcode",
    group_col: str = "source",
    thresholds: Optional[str] = None,
    max_level: int = 17,
) -> pd.DataFrame:
    inferred = infer_max_level_from_lincode(df[lin_col])
    max_level = min(int(max_level), int(inferred))
    ks = parse_thresholds(thresholds, max_level)

    rows = []
    for k in ks:
        prefixes = df[lin_col].astype(str).str.split("_", expand=True).iloc[:, :k]
        prefix_code = prefixes.apply(lambda r: "_".join(r.values.astype(str)), axis=1)

        if group_col:
            for grp, sub in df.groupby(group_col, dropna=False):
                sub_prefix = prefix_code.loc[sub.index]
                rows.append({"threshold": k, group_col: grp, "n_clusters": sub_prefix.nunique()})
        else:
            rows.append({"threshold": k, "n_clusters": prefix_code.nunique()})

    out = pd.DataFrame(rows)
    return out


def run_diversification(
    inpath: str,
    outdir: str,
    lin_col: str = "LINcode",
    group_col: str = "source",
    thresholds: Optional[str] = None,
    max_level: int = 17,
    title: Optional[str] = None,
    formats: str = "png,svg",
) -> None:
    outdir_p = Path(outdir)
    outdir_p.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(inpath, sep="\t", low_memory=False)
    div = lin_diversification(
        df,
        lin_col=lin_col,
        group_col=group_col,
        thresholds=thresholds,
        max_level=max_level,
    )

    div.to_csv(outdir_p / "diversification.tsv", sep="\t", index=False)

    max_level = int(min(max_level, div["threshold"].max())) if not div.empty else int(max_level)
    fmts = [f.strip() for f in formats.split(',') if f.strip()]

    # Reservoir-focused view (exclude human + other by default)
    if group_col in div.columns:
        res = div[~div[group_col].astype(str).str.lower().isin(["human", "other"])].copy()
        plot_diversification_curve(
            res,
            group_col=group_col,
            max_level=max_level,
            title=title or "Diversification by LIN threshold (reservoir sources)",
            outpath=outdir_p / "diversification",
            formats=fmts,
        )

        all_ = div.copy()
        plot_diversification_curve(
            all_,
            group_col=group_col,
            max_level=max_level,
            title=title or "Diversification by LIN threshold (all sources)",
            outpath=outdir_p / "diversification_all_sources",
            formats=fmts,
        )
    else:
        plot_diversification_curve(
            div,
            group_col=None,
            max_level=max_level,
            title=title or "Diversification by LIN threshold",
            outpath=outdir_p / "diversification",
            formats=fmts,
        )

