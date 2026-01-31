# linwalker/plotting.py

from __future__ import annotations
import matplotlib.pyplot as plt
import pandas as pd
from .palette import SOURCE_COLOURS

# Stable ordering for attribution-style plots
SOURCE_ORDER = ["chicken", "pig", "ruminant", "wild bird", "human", "other"]


def _sorted_sources(sources):
    srcs = [str(s).lower() for s in sources]
    order = {k: i for i, k in enumerate(SOURCE_ORDER)}
    return sorted(srcs, key=lambda x: order.get(x, 999))


def plot_diversification(
    div_df: pd.DataFrame,
    title: str | None = None,
    outpath: str | None = None,
    exclude_sources: set[str] | None = None,
):
    """
    Publication-style diversification figure (unique LIN IDs vs threshold by source).

    Note
    ----
    Human/other often dominate counts and can compress reservoir curves.
    By default we exclude {"human","other"} unless override provided.
    """
    if exclude_sources is None:
        exclude_sources = {"human", "other"}

    fig = plt.figure(figsize=(10, 7))
    for src in _sorted_sources(div_df["group"].unique()):
        if src in exclude_sources:
            continue
        sub = div_df[div_df["group"].astype(str).str.lower() == src]
        if sub.empty:
            continue
        plt.plot(
            sub["LIN_level"].astype(int),
            sub["n_unique_LINs"],
            label=src,
            color=SOURCE_COLOURS.get(src, "#000000"),
            linewidth=3.2,
        )

    plt.xlabel("LIN threshold (1–17)")
    plt.ylabel("Number of unique LIN IDs")
    if title:
        plt.title(title)
    plt.legend(title="Source", frameon=True)
    plt.grid(True, linestyle="--", alpha=0.35)

    # clean ticks: show every 2 levels
    plt.xticks(list(range(1, 18, 2)))
    plt.tight_layout()
    if outpath:
        plt.savefig(outpath, dpi=300)
    return fig


def plot_mixed_species(mix_df: pd.DataFrame, title: str | None = None, outpath: str | None = None):
    """
    Publication-style introgression curve: proportion of mixed-species LIN clusters vs threshold.
    """
    fig = plt.figure(figsize=(10, 7))
    plt.plot(mix_df["LIN_level"].astype(int), mix_df["prop_mixed_species"], linewidth=3.2)
    plt.xlabel("LIN threshold (1–17)")
    plt.ylabel("Proportion of mixed-species LIN clusters")
    if title:
        plt.title(title)
    plt.grid(True, linestyle="--", alpha=0.35)
    plt.xticks(list(range(1, 18, 2)))
    plt.ylim(0, min(1.0, max(0.05, float(mix_df["prop_mixed_species"].max()) * 1.05)))
    plt.tight_layout()
    if outpath:
        plt.savefig(outpath, dpi=300)
    return fig


def plot_lsdd_by_source(df_lsdd: pd.DataFrame, source_col: str = "source", title: str | None = None, outpath: str | None = None):
    """
    Boxplot of LSDD stratified by source, coloured using SOURCE_COLOURS.
    """
    fig = plt.figure(figsize=(10, 7))
    df = df_lsdd.copy()
    df[source_col] = df[source_col].astype(str).str.lower()

    sources = _sorted_sources(df[source_col].unique())
    data = [pd.to_numeric(df.loc[df[source_col] == s, "LSDD"], errors="coerce").dropna().values for s in sources]

    bp = plt.boxplot(data, labels=sources, vert=True, showfliers=False, patch_artist=True)

    for patch, s in zip(bp["boxes"], sources):
        patch.set_facecolor(SOURCE_COLOURS.get(s, "#CCCCCC"))
        patch.set_edgecolor("#000000")
        patch.set_linewidth(1.4)

    for element in ["whiskers", "caps", "medians"]:
        for line in bp[element]:
            line.set_color("#000000")
            line.set_linewidth(1.4)

    plt.ylabel("LSDD (earliest discordant LIN level)")
    if title:
        plt.title(title)
    plt.grid(True, linestyle="--", alpha=0.35, axis="y")
    plt.xticks(rotation=20, ha="right")
    plt.tight_layout()
    if outpath:
        plt.savefig(outpath, dpi=300)
    return fig



def plot_stcc_concordance(stcc_df: pd.DataFrame, title: str | None = None, outpath: str | None = None):
    """
    Plot LIN threshold vs concordance metrics for MLST ST and clonal complex.
    Shows proportion of pure clusters and weighted purity for both ST and CC.
    """
    fig = plt.figure(figsize=(10, 7))

    x = stcc_df["LIN_level"].astype(int)

    # Proportion pure
    plt.plot(x, stcc_df["prop_pure_ST"], label="ST purity (prop pure)", linewidth=3.2)
    plt.plot(x, stcc_df["prop_pure_CC"], label="CC purity (prop pure)", linewidth=3.2, linestyle="--")

    # Weighted purity (secondary axis)
    ax1 = plt.gca()
    ax2 = ax1.twinx()
    ax2.plot(x, stcc_df["weighted_purity_ST"], label="ST purity (weighted)", linewidth=3.2, alpha=0.85)
    ax2.plot(x, stcc_df["weighted_purity_CC"], label="CC purity (weighted)", linewidth=3.2, alpha=0.85, linestyle="--")

    ax1.set_xlabel("LIN threshold (1–17)")
    ax1.set_ylabel("Proportion of pure clusters")
    ax2.set_ylabel("Weighted purity (max fraction within cluster)")

    if title:
        plt.title(title)

    ax1.grid(True, linestyle="--", alpha=0.35)
    ax1.set_xticks(list(range(1, 18, 2)))
    ax1.set_ylim(0, 1.0)
    ax2.set_ylim(0, 1.0)

    # Combined legend
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, frameon=True, loc="lower right")

    plt.tight_layout()
    if outpath:
        plt.savefig(outpath, dpi=300)
    return fig
