"""Outbreak / public-health summaries for LINwalker.

This module is intentionally *descriptive* (no inference). It produces plots and
TSVs aimed at outbreak-style questions:

- How do per-isolate LIN cluster sizes change with LIN threshold?
- Which clusters dominate at a chosen threshold?
- Do the biggest clusters mix sources/species/countries?
- Optional epi-curve by date (if a date column is present)

Inputs are typically the `*_LINwalker_min.tsv` table produced by `linwalker prep`,
but any TSV with a full LIN code column will work.
"""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np
import pandas as pd

from .plotting import save_figure


def _ensure_outdir(outdir: str | Path) -> Path:
    p = Path(outdir)
    p.mkdir(parents=True, exist_ok=True)
    return p


def lin_prefix(lincode: str, level: int, sep: str = "_") -> str:
    """Return the LIN prefix at 1-indexed `level`.

    LIN codes may be stored with '_' separators (typical) or '.' or other.
    We standardise by splitting on underscores first; if not present, try '.';
    else treat as already tokenised by whitespace.
    """
    if pd.isna(lincode):
        return ""
    s = str(lincode)
    if "_" in s:
        parts = s.split("_")
    elif "." in s:
        parts = s.split(".")
    else:
        parts = s.split()
    level = max(1, min(level, len(parts)))
    return sep.join(parts[:level])


def infer_max_level(series: pd.Series) -> int:
    """Infer maximum LIN depth from a series of LIN codes."""
    s = series.dropna().astype(str)
    if s.empty:
        return 0
    # assume underscore separated
    counts = s.str.count("_") + 1
    return int(counts.max())


@dataclass
class ClusterSummary:
    level: int
    n_isolates: int
    n_clusters: int
    median_cluster_size: float
    p25: float
    p75: float
    max_cluster_size: int


def per_isolate_cluster_sizes(
    df: pd.DataFrame,
    lin_col: str = "LINcode",
    levels: Optional[Iterable[int]] = None,
) -> Dict[int, np.ndarray]:
    """Return per-isolate cluster sizes for each LIN level.

    For each isolate, we compute the size of its LIN-prefix cluster at that level.
    This yields a distribution of sizes (length n isolates) per level.
    """
    if lin_col not in df.columns:
        raise ValueError(f"Missing LIN column: {lin_col}")

    max_level = infer_max_level(df[lin_col])
    if max_level == 0:
        raise ValueError("Could not infer LIN depth (no LIN codes?)")

    if levels is None:
        levels = range(1, max_level + 1)
    else:
        levels = [int(x) for x in levels]

    out: Dict[int, np.ndarray] = {}

    lins = df[lin_col].astype(str)
    for lvl in levels:
        prefixes = lins.map(lambda x: lin_prefix(x, lvl))
        sizes = prefixes.map(prefixes.value_counts()).to_numpy()
        out[lvl] = sizes

    return out


def summarise_levels(size_map: Dict[int, np.ndarray]) -> List[ClusterSummary]:
    rows: List[ClusterSummary] = []
    for lvl, sizes in sorted(size_map.items(), key=lambda x: x[0]):
        if len(sizes) == 0:
            continue
        rows.append(
            ClusterSummary(
                level=int(lvl),
                n_isolates=int(len(sizes)),
                n_clusters=int(len(np.unique(sizes)) if False else 0),  # placeholder, filled below
                median_cluster_size=float(np.median(sizes)),
                p25=float(np.percentile(sizes, 25)),
                p75=float(np.percentile(sizes, 75)),
                max_cluster_size=int(np.max(sizes)),
            )
        )

    # n_clusters can't be derived from sizes alone reliably; leave 0 here.
    return rows


def top_clusters_table(
    df: pd.DataFrame,
    level: int,
    lin_col: str = "LINcode",
    id_col: str = "id",
    group_cols: Optional[List[str]] = None,
    top_n: int = 25,
) -> pd.DataFrame:
    """Return a table of top LIN clusters at a given level.

    Includes size and (optionally) composition summaries for given columns.
    """
    if group_cols is None:
        group_cols = [c for c in ["species", "source", "country"] if c in df.columns]

    tmp = df.copy()
    tmp["lin_prefix"] = tmp[lin_col].astype(str).map(lambda x: lin_prefix(x, level))

    counts = tmp["lin_prefix"].value_counts().rename("n").reset_index()
    counts = counts.rename(columns={"index": "lin_prefix"}).head(top_n)

    # composition: for each group col, add the top label + proportion
    out = counts.copy()
    for col in group_cols:
        comp = (
            tmp.groupby("lin_prefix")[col]
            .apply(lambda s: s.value_counts(dropna=False).head(3))
            .rename("count")
            .reset_index()
        )
        # turn into compact "label:count" string per prefix
        comp["label"] = comp[col].astype(str) + ":" + comp["count"].astype(str)
        comp2 = comp.groupby("lin_prefix")["label"].apply(lambda x: ";".join(x)).reset_index()
        out = out.merge(comp2, on="lin_prefix", how="left")
        out = out.rename(columns={"label": f"top_{col}"})

    return out


def parse_date_safe(x: str) -> Optional[datetime]:
    if pd.isna(x):
        return None
    s = str(x).strip()
    if not s:
        return None
    # accept YYYY-MM-DD or YYYY/MM/DD or YYYY
    for fmt in ("%Y-%m-%d", "%Y/%m/%d", "%Y-%m", "%Y/%m", "%Y"):
        try:
            return datetime.strptime(s, fmt)
        except ValueError:
            continue
    return None


def epi_curve(
    df: pd.DataFrame,
    date_col: str,
    outdir: str | Path,
    title: str = "Epidemic curve",
    formats: Tuple[str, ...] = ("png", "svg"),
):
    """Simple epi-curve plot (counts by month) if dates are parseable."""
    import matplotlib.pyplot as plt

    outdir = _ensure_outdir(outdir)
    if date_col not in df.columns:
        raise ValueError(f"Missing date column: {date_col}")

    dt = df[date_col].map(parse_date_safe)
    dt = dt.dropna()
    if dt.empty:
        raise ValueError("No parseable dates found")

    # bin by month
    months = dt.map(lambda d: d.strftime("%Y-%m"))
    counts = months.value_counts().sort_index()

    fig, ax = plt.subplots(figsize=(12, 4.5))
    ax.bar(counts.index, counts.values)
    ax.set_title(title)
    ax.set_ylabel("Isolates")
    ax.set_xlabel("Month")
    ax.tick_params(axis="x", labelrotation=45)
    fig.tight_layout()

    save_figure(fig, outdir / "epi_curve", formats=formats)
    plt.close(fig)


def plot_cluster_size_summary(
    size_map: Dict[int, np.ndarray],
    outdir: str | Path,
    title: str = "Per-isolate LIN cluster sizes across thresholds",
    formats: Tuple[str, ...] = ("png", "svg"),
    ylog: bool = True,
):
    """Plot median + IQR of per-isolate cluster size as threshold increases."""
    import matplotlib.pyplot as plt

    outdir = _ensure_outdir(outdir)

    levels = np.array(sorted(size_map.keys()), dtype=int)
    med = np.array([np.median(size_map[l]) for l in levels], dtype=float)
    p25 = np.array([np.percentile(size_map[l], 25) for l in levels], dtype=float)
    p75 = np.array([np.percentile(size_map[l], 75) for l in levels], dtype=float)

    fig, ax = plt.subplots(figsize=(12, 4.5))
    ax.plot(levels, med, marker="o")
    ax.fill_between(levels, p25, p75, alpha=0.2)
    ax.set_title(title)
    ax.set_xlabel("LIN threshold")
    ax.set_ylabel("Cluster size (per isolate)")
    ax.set_xticks(levels)
    if ylog:
        ax.set_yscale("log")
    fig.tight_layout()

    save_figure(fig, outdir / "cluster_size_summary", formats=formats)
    plt.close(fig)


def plot_cluster_size_boxplot(
    size_map: Dict[int, np.ndarray],
    outdir: str | Path,
    title: str = "Per-isolate cluster size distribution by LIN threshold",
    formats: Tuple[str, ...] = ("png", "svg"),
    show_points: bool = True,
    max_points: int = 20000,
    ylog: bool = True,
):
    """Boxplots across thresholds, with optional downsampled per-isolate points."""
    import matplotlib.pyplot as plt

    outdir = _ensure_outdir(outdir)

    levels = sorted(size_map.keys())
    data = [size_map[l] for l in levels]

    fig, ax = plt.subplots(figsize=(13, 5))
    ax.boxplot(data, labels=[str(l) for l in levels], showfliers=False)

    if show_points:
        rng = np.random.default_rng(1)
        for i, arr in enumerate(data, start=1):
            arr = np.asarray(arr)
            if arr.size > max_points:
                idx = rng.choice(arr.size, size=max_points, replace=False)
                arr = arr[idx]
            x = rng.normal(loc=i, scale=0.07, size=arr.size)
            ax.scatter(x, arr, s=2, alpha=0.15)

    ax.set_title(title)
    ax.set_xlabel("LIN threshold")
    ax.set_ylabel("Cluster size (per isolate)")
    if ylog:
        ax.set_yscale("log")
    fig.tight_layout()

    save_figure(fig, outdir / "cluster_size_boxplot", formats=formats)
    plt.close(fig)

