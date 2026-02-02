from __future__ import annotations

from pathlib import Path
from typing import Iterable, List, Optional

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

from .palette import get_palette, RESERVOIR_SOURCES


def _resolve_outbase(
    *,
    outpath_base: Optional[Path] = None,
    outdir: Optional[Path] = None,
    filename: Optional[str] = None,
) -> Path:
    if outpath_base is not None:
        return Path(outpath_base)
    if outdir is None or filename is None:
        raise ValueError("Provide either outpath_base or (outdir + filename)")
    return Path(outdir) / filename


def _save_fig(fig: plt.Figure, outpath_base: Path, formats: Iterable[str], dpi: int = 300) -> None:
    outpath_base = Path(outpath_base)
    outpath_base.parent.mkdir(parents=True, exist_ok=True)
    for fmt in formats:
        fig.savefig(str(outpath_base.with_suffix(f".{fmt}")), dpi=dpi, bbox_inches="tight")


def plot_diversification_curve(
    div_table: pd.DataFrame,
    *,
    outdir: Optional[Path] = None,
    filename: str = "diversification",
    outpath_base: Optional[Path] = None,
    formats: Optional[List[str]] = None,
    title: str = "Unique LIN IDs vs. LIN threshold by source",
    max_level: Optional[int] = None,
    include_sources: Optional[List[str]] = None,
) -> None:
    """Plot number of unique LIN prefixes vs LIN level, stratified by source."""
    if formats is None:
        formats = ["png", "svg"]
    outbase = _resolve_outbase(outpath_base=outpath_base, outdir=outdir, filename=filename)

    df = div_table.copy()
    if "lin_level" not in df.columns:
        # Back-compat: accept "threshold" or "LIN_level"
        if "threshold" in df.columns:
            df = df.rename(columns={"threshold": "lin_level"})
        elif "LIN_level" in df.columns:
            df = df.rename(columns={"LIN_level": "lin_level"})
    df["lin_level"] = pd.to_numeric(df["lin_level"], errors="coerce")
    df = df.dropna(subset=["lin_level"])

    if max_level is None:
        max_level = int(df["lin_level"].max()) if len(df) else 17

    if include_sources is not None:
        df = df[df["source"].isin(include_sources)].copy()

    palette = get_palette()

    fig = plt.figure(figsize=(11, 7))
    ax = fig.add_subplot(111)

    for src, sub in df.groupby("source"):
        sub = sub.sort_values("lin_level")
        ax.plot(sub["lin_level"], sub["n_unique"], linewidth=3, label=src, color=palette.get(src, None))

    ax.set_title(title)
    ax.set_xlabel(f"LIN threshold (1–{max_level})")
    ax.set_ylabel("Number of unique LIN IDs")
    ax.set_xlim(1, max_level)
    # Avoid "squashed" tick labels: show only a few major labels.
    # User-requested majors: 1, 5, 10, 15.
    major_ticks = [t for t in (1, 5, 10, 15) if 1 <= t <= max_level]
    ax.set_xticks(major_ticks)
    # Keep minor ticks at every level for vertical guide lines.
    ax.set_xticks(list(range(1, max_level + 1)), minor=True)
    ax.grid(True, which="major", linestyle="--", alpha=0.3)
    ax.grid(True, which="minor", axis="x", linestyle="--", alpha=0.12)
    for label in ax.get_xticklabels():
        label.set_rotation(0)
        label.set_horizontalalignment("center")

    ax.legend(title="Source", frameon=True)

    _save_fig(fig, outbase, formats=formats)
    plt.close(fig)


def plot_mixed_species(
    mixed_table: pd.DataFrame,
    *,
    outdir: Optional[Path] = None,
    filename: str = "mixed_species",
    outpath_base: Optional[Path] = None,
    formats: Optional[List[str]] = None,
    title: str = "Mixed-species LIN prefixes vs LIN threshold",
    max_level: Optional[int] = None,
) -> None:
    if formats is None:
        formats = ["png", "svg"]
    outbase = _resolve_outbase(outpath_base=outpath_base, outdir=outdir, filename=filename)

    df = mixed_table.copy()
    if "lin_level" not in df.columns and "threshold" in df.columns:
        df = df.rename(columns={"threshold": "lin_level"})
    df["lin_level"] = pd.to_numeric(df["lin_level"], errors="coerce")
    df = df.dropna(subset=["lin_level"])

    if max_level is None:
        max_level = int(df["lin_level"].max()) if len(df) else 17

    fig = plt.figure(figsize=(11, 7))
    ax = fig.add_subplot(111)

    df = df.sort_values("lin_level")
    ax.plot(df["lin_level"], df["mixed_fraction"], linewidth=3)

    ax.set_title(title)
    ax.set_xlabel(f"LIN threshold (1–{max_level})")
    ax.set_ylabel("Fraction of prefixes shared by >1 species")
    ax.set_xlim(1, max_level)
    ax.set_xticks(list(range(1, max_level + 1)))
    ax.grid(True, linestyle="--", alpha=0.3)

    _save_fig(fig, outbase, formats=formats)
    plt.close(fig)


def plot_lsdd(
    lsdd_table: pd.DataFrame,
    *,
    outdir: Optional[Path] = None,
    filename: str = "lsdd",
    outpath_base: Optional[Path] = None,
    formats: Optional[List[str]] = None,
    title: str = "LSDD (within vs between species) by LIN threshold",
    max_level: Optional[int] = None,
) -> None:
    if formats is None:
        formats = ["png", "svg"]
    outbase = _resolve_outbase(outpath_base=outpath_base, outdir=outdir, filename=filename)

    df = lsdd_table.copy()
    if "lin_level" not in df.columns and "threshold" in df.columns:
        df = df.rename(columns={"threshold": "lin_level"})
    df["lin_level"] = pd.to_numeric(df["lin_level"], errors="coerce")
    df = df.dropna(subset=["lin_level"])

    if max_level is None:
        max_level = int(df["lin_level"].max()) if len(df) else 17

    fig = plt.figure(figsize=(11, 7))
    ax = fig.add_subplot(111)

    df = df.sort_values("lin_level")
    ax.plot(df["lin_level"], df["lsdd"], linewidth=3)

    ax.set_title(title)
    ax.set_xlabel(f"LIN threshold (1–{max_level})")
    ax.set_ylabel("LSDD")
    ax.set_xlim(1, max_level)
    ax.set_xticks(list(range(1, max_level + 1)))
    ax.grid(True, linestyle="--", alpha=0.3)

    _save_fig(fig, outbase, formats=formats)
    plt.close(fig)


def plot_stcc_concordance(
    stcc_table: pd.DataFrame,
    *,
    outdir: Optional[Path] = None,
    filename: str = "stcc_concordance",
    outpath_base: Optional[Path] = None,
    formats: Optional[List[str]] = None,
    title: str = "ST/CC concordance vs LIN threshold",
    max_level: Optional[int] = None,
) -> None:
    if formats is None:
        formats = ["png", "svg"]
    outbase = _resolve_outbase(outpath_base=outpath_base, outdir=outdir, filename=filename)

    df = stcc_table.copy()
    if "lin_level" not in df.columns and "threshold" in df.columns:
        df = df.rename(columns={"threshold": "lin_level"})
    df["lin_level"] = pd.to_numeric(df["lin_level"], errors="coerce")
    df = df.dropna(subset=["lin_level"])

    if max_level is None:
        max_level = int(df["lin_level"].max()) if len(df) else 17

    fig = plt.figure(figsize=(11, 7))
    ax = fig.add_subplot(111)

    df = df.sort_values("lin_level")
    st_col = "st_purity" if "st_purity" in df.columns else "mean_purity_ST"
    cc_col = "cc_purity" if "cc_purity" in df.columns else ("mean_purity_CC" if "mean_purity_CC" in df.columns else None)

    # Matplotlib can't cast pandas.NA (NAType) to float; coerce to numeric so
    # missing values become NaN.
    if st_col in df.columns:
        df[st_col] = pd.to_numeric(df[st_col], errors="coerce")
    if cc_col and cc_col in df.columns:
        df[cc_col] = pd.to_numeric(df[cc_col], errors="coerce")

    plotted_any = False

    if st_col in df.columns:
        d1 = df[["lin_level", st_col]].dropna()
        if len(d1):
            ax.plot(d1["lin_level"], d1[st_col], linewidth=3, label="ST purity")
            plotted_any = True

    if cc_col is not None and cc_col in df.columns:
        d2 = df[["lin_level", cc_col]].dropna()
        if len(d2):
            ax.plot(d2["lin_level"], d2[cc_col], linewidth=3, label="CC purity")
            plotted_any = True

    if not plotted_any:
        ax.text(
            0.5,
            0.5,
            "No ST/CC concordance metrics available\n(check ST/CC columns or run prep)",
            ha="center",
            va="center",
            transform=ax.transAxes,
            fontsize=14,
        )

    ax.set_title(title)
    ax.set_xlabel(f"LIN threshold (1–{max_level})")
    ax.set_ylabel("Purity (1.0 = perfect concordance)")
    ax.set_xlim(1, max_level)
    ax.set_xticks(list(range(1, max_level + 1)))
    ax.set_ylim(0, 1.05)
    ax.grid(True, linestyle="--", alpha=0.3)
    ax.legend(frameon=True)

    _save_fig(fig, outbase, formats=formats)
    plt.close(fig)


def plot_outbreak_top_clusters(
    top_clusters: pd.DataFrame,
    *,
    outdir: Optional[Path] = None,
    filename: str = "top_clusters",
    outpath_base: Optional[Path] = None,
    formats: Optional[List[str]] = None,
    title: str = "Top LIN clusters",
) -> None:
    if formats is None:
        formats = ["png", "svg"]
    outbase = _resolve_outbase(outpath_base=outpath_base, outdir=outdir, filename=filename)

    df = top_clusters.copy()
    if df.empty:
        return

    fig = plt.figure(figsize=(11, 7))
    ax = fig.add_subplot(111)

    df = df.sort_values("n", ascending=True)
    col = "lin_prefix" if "lin_prefix" in df.columns else "lin_cluster"
    ax.barh(df[col].astype(str), df["n"].astype(int))

    ax.set_title(title)
    ax.set_xlabel("Isolate count")
    ax.set_ylabel("LIN prefix")
    ax.grid(True, axis="x", linestyle="--", alpha=0.3)

    _save_fig(fig, outbase, formats=formats)
    plt.close(fig)


def plot_outbreak_epicurve(
    date_counts: pd.DataFrame,
    *,
    outdir: Optional[Path] = None,
    filename: str = "epicurve",
    outpath_base: Optional[Path] = None,
    formats: Optional[List[str]] = None,
    title: str = "Counts over time",
) -> None:
    if formats is None:
        formats = ["png", "svg"]
    outbase = _resolve_outbase(outpath_base=outpath_base, outdir=outdir, filename=filename)

    df = date_counts.copy()
    if df.empty or "date" not in df.columns:
        return

    df["date"] = pd.to_datetime(df["date"], errors="coerce")
    df = df.dropna(subset=["date"]).sort_values("date")

    fig = plt.figure(figsize=(11, 7))
    ax = fig.add_subplot(111)

    ax.plot(df["date"], df["n"], linewidth=3)

    ax.set_title(title)
    ax.set_xlabel("Date")
    ax.set_ylabel("Isolate count")
    ax.grid(True, linestyle="--", alpha=0.3)

    _save_fig(fig, outbase, formats=formats)
    plt.close(fig)
