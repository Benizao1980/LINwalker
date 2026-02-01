from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import pandas as pd

from .utils import lin_prefix


@dataclass
class OutbreakResult:
    top_clusters: pd.DataFrame
    cluster_source_counts: pd.DataFrame
    cluster_country_counts: Optional[pd.DataFrame]
    cluster_date_counts: Optional[pd.DataFrame]


def outbreak_descriptives(
    df: pd.DataFrame,
    *,
    lin_col: str = "lin_code",
    source_col: str = "source",
    country_col: Optional[str] = "country",
    date_col: Optional[str] = None,
    top_threshold: int = 12,
    top_n: int = 25,
) -> OutbreakResult:
    """Descriptive outbreak/public-health summaries from LIN codes.

    Creates:
    - top LIN clusters at a given threshold (counts)
    - source composition for those clusters
    - optional country and date summaries if columns are available
    """

    df = df.copy()
    df["lin_cluster"] = df[lin_col].astype(str).map(lambda x: lin_prefix(x, top_threshold))

    top = (
        df["lin_cluster"].value_counts(dropna=False)
        .rename_axis("lin_cluster")
        .reset_index(name="n")
        .head(top_n)
    )

    # Source composition for top clusters
    top_set = set(top["lin_cluster"].tolist())
    dft = df[df["lin_cluster"].isin(top_set)].copy()

    src_counts = (
        dft.groupby(["lin_cluster", source_col], dropna=False)
        .size()
        .reset_index(name="n")
        .sort_values(["lin_cluster", "n"], ascending=[True, False])
    )

    country_counts = None
    if country_col and country_col in df.columns:
        country_counts = (
            dft.groupby(["lin_cluster", country_col], dropna=False)
            .size()
            .reset_index(name="n")
            .sort_values(["lin_cluster", "n"], ascending=[True, False])
        )

    date_counts = None
    if date_col and date_col in df.columns:
        # attempt to parse; keep YYYY-MM if possible
        dates = pd.to_datetime(dft[date_col], errors="coerce")
        dft = dft.assign(_ym=dates.dt.to_period("M").astype(str))
        date_counts = (
            dft.groupby(["lin_cluster", "_ym"], dropna=False)
            .size()
            .reset_index(name="n")
            .rename(columns={"_ym": "year_month"})
            .sort_values(["lin_cluster", "year_month"])
        )

    return OutbreakResult(
        top_clusters=top,
        cluster_source_counts=src_counts,
        cluster_country_counts=country_counts,
        cluster_date_counts=date_counts,
    )
