"""Plotting utilities.

Design goals:
- deterministic axes (LIN thresholds are integer levels)
- deterministic source ordering/colours
- light dependency footprint (matplotlib only)

All plotting functions save both PNG and SVG by default.
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable, List, Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .palette import SOURCE_COLOURS

# Stable ordering (and a sane ecological default)
SOURCE_ORDER = ["chicken", "ruminant", "pig", "wild bird", "human", "other"]


def _sorted_sources(sources: Iterable[str]) -> List[str]:
    srcs = [str(s).lower() for s in sources]
    order = {k: i for i, k in enumerate(SOURCE_ORDER)}
    return sorted(set(srcs), key=lambda s: order.get(s, 999))


def _ensure_outdir(outdir: str | Path) -> Path:
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    return outdir


def _configure_axis(ax, max_level: int, title: Optional[str], ylabel: str):
    ax.set_xlabel("LIN threshold")
    ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title)

    # LIN thresholds are discrete 1..max_level
    ax.set_xlim(1, max_level)
    ax.set_xticks(list(range(1, max_level + 1)))
    ax.grid(True, alpha=0.25)


def _save(fig, outbase: Path, formats: Iterable[str] = ("png", "svg")):
    for fmt in formats:
        fmt = fmt.lower().lstrip(".")
        fig.savefig(outbase.with_suffix(f".{fmt}"), bbox_inches="tight", dpi=300)


def save_figure(fig, outbase: Path, formats: Iterable[str] = ("png", "svg")):
    """Public wrapper used by optional modules (e.g. outbreak)."""
    _save(fig, outbase, formats=formats)


def plot_diversification_curve(
    df: pd.DataFrame,
    outdir: str | Path,
    title: str = "Diversification by LIN threshold",
    max_level: Optional[int] = None,
    formats: Iterable[str] = ("png", "svg"),
    filename: str = "diversification",
):
    """Line plot of unique clusters vs LIN threshold, stratified by source."""

    outdir = _ensure_outdir(outdir)
    d = df.copy()
    d["threshold"] = pd.to_numeric(d["threshold"], errors="coerce").astype("Int64")
    d = d.dropna(subset=["threshold"])

    if max_level is None:
        max_level = int(d["threshold"].max())

    fig, ax = plt.subplots(figsize=(10, 4))
    for src in _sorted_sources(d["source"].unique()):
        sub = d[d["source"].str.lower() == src].sort_values("threshold")
        if sub.empty:
            continue
        colour = SOURCE_COLOURS.get(src, "#777777")
        ax.plot(sub["threshold"].to_numpy(), sub["n_unique"].to_numpy(), label=src, linewidth=2)

    _configure_axis(ax, max_level=max_level, title=title, ylabel="Unique LIN clusters")
    ax.legend(frameon=False, ncol=3)

    _save(fig, outdir / filename, formats=formats)
    plt.close(fig)


def plot_mixed_species(
    df: pd.DataFrame,
    outdir: str | Path,
    title: str = "Mixed-species LIN clusters by threshold",
    max_level: Optional[int] = None,
    formats: Iterable[str] = ("png", "svg"),
    filename: str = "mixed_species",
):
    """Proportion of LIN clusters that contain >1 species, by LIN threshold."""

    outdir = _ensure_outdir(outdir)
    d = df.copy()
    d["threshold"] = pd.to_numeric(d["threshold"], errors="coerce").astype("Int64")
    d = d.dropna(subset=["threshold"])

    if max_level is None:
        max_level = int(d["threshold"].max())

    fig, ax = plt.subplots(figsize=(10, 4))
    ax.plot(d["threshold"].to_numpy(), d["prop_mixed"].to_numpy(), linewidth=2)

    _configure_axis(ax, max_level=max_level, title=title, ylabel="Proportion mixed-species")
    ax.set_ylim(0, 1)

    _save(fig, outdir / filename, formats=formats)
    plt.close(fig)


def plot_lsdd_by_source(
    df: pd.DataFrame,
    outdir: str | Path,
    title: str = "LSDD by source",
    max_level: int = 17,
    formats: Iterable[str] = ("png", "svg"),
    filename: str = "lsdd_by_source",
    show_points: bool = True,
):
    """Boxplot of LSDD distributions by source, with optional jittered points."""

    outdir = _ensure_outdir(outdir)

    d = df.copy()
    d["source"] = d["source"].astype(str).str.lower()
    d["lsdd"] = pd.to_numeric(d["lsdd"], errors="coerce")
    d = d.dropna(subset=["lsdd"])

    sources = _sorted_sources(d["source"].unique())
    data = [d.loc[d["source"] == s, "lsdd"].to_numpy() for s in sources]

    fig, ax = plt.subplots(figsize=(10, 4))
    bp = ax.boxplot(
        data,
        labels=sources,
        showfliers=False,
        patch_artist=True,
        widths=0.65,
    )

    # colour boxes by source
    for patch, s in zip(bp["boxes"], sources):
        patch.set_facecolor(SOURCE_COLOURS.get(s, "#DDDDDD"))
        patch.set_alpha(0.6)

    if show_points:
        for i, s in enumerate(sources, start=1):
            y = d.loc[d["source"] == s, "lsdd"].to_numpy()
            if len(y) == 0:
                continue
            x = np.random.normal(i, 0.06, size=len(y))
            ax.scatter(x, y, s=8, alpha=0.35, linewidths=0)

    ax.set_ylabel("LSDD (earliest discordant LIN level)")
    ax.set_ylim(0.5, max_level + 0.5)
    ax.set_yticks(list(range(1, max_level + 1)))
    if title:
        ax.set_title(title)
    ax.grid(True, axis="y", alpha=0.25)

    _save(fig, outdir / filename, formats=formats)
    plt.close(fig)


def plot_stcc_concordance(
    df: pd.DataFrame,
    outdir: str | Path,
    title: str = "LIN ↔ MLST concordance by threshold",
    max_level: Optional[int] = None,
    formats: Iterable[str] = ("png", "svg"),
    filename: str = "stcc_concordance",
):
    """Plot ST/CC purity metrics vs LIN threshold."""

    outdir = _ensure_outdir(outdir)
    d = df.copy()
    d["threshold"] = pd.to_numeric(d["threshold"], errors="coerce").astype("Int64")
    d = d.dropna(subset=["threshold"])

    if max_level is None:
        max_level = int(d["threshold"].max())

    fig, ax = plt.subplots(figsize=(10, 4))

    for col, label in [
        ("prop_pure_st", "Pure-ST clusters"),
        ("weighted_purity_st", "Weighted ST purity"),
        ("prop_pure_cc", "Pure-CC clusters"),
        ("weighted_purity_cc", "Weighted CC purity"),
    ]:
        if col in d.columns:
            ax.plot(d["threshold"].to_numpy(), d[col].to_numpy(), label=label, linewidth=2)

    _configure_axis(ax, max_level=max_level, title=title, ylabel="Concordance")
    ax.set_ylim(0, 1)
    ax.legend(frameon=False, ncol=2)

    _save(fig, outdir / filename, formats=formats)
    plt.close(fig)
