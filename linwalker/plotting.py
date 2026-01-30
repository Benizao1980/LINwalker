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


def plot_diversification(div_df: pd.DataFrame, title: str | None = None, outpath: str | None = None):
    """
    Publication-style diversification figure (unique LIN IDs vs threshold by source).
    """
    fig = plt.figure(figsize=(10, 7))
    for src in _sorted_sources(div_df["group"].unique()):
        sub = div_df[div_df["group"].astype(str).str.lower() == src]
        if sub.empty:
            continue
        plt.plot(
            sub["LIN_level"],
            sub["n_unique_LINs"],
            label=src,
            color=SOURCE_COLOURS.get(src, "#000000"),
            linewidth=3.2
        )

    plt.xlabel("LIN threshold (1–17)")
    plt.ylabel("Number of unique LIN IDs")
    if title:
        plt.title(title)
    plt.legend(title="Source", frameon=True)
    plt.grid(True, linestyle="--", alpha=0.35)
    plt.xticks(range(int(div_df["LIN_level"].min()), int(div_df["LIN_level"].max()) + 1, 2))
    plt.tight_layout()
    if outpath:
        plt.savefig(outpath, dpi=300)
    return fig


def plot_mixed_species(mix_df: pd.DataFrame, title: str | None = None, outpath: str | None = None):
    """
    Publication-style introgression curve: proportion of mixed-species LIN clusters vs threshold.
    """
    fig = plt.figure(figsize=(10, 7))
    plt.plot(mix_df["LIN_level"], mix_df["prop_mixed_species"], linewidth=3.2)
    plt.xlabel("LIN threshold (1–17)")
    plt.ylabel("Proportion of mixed-species LIN clusters")
    if title:
        plt.title(title)
    plt.grid(True, linestyle="--", alpha=0.35)
    plt.xticks(range(int(mix_df["LIN_level"].min()), int(mix_df["LIN_level"].max()) + 1, 2))
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
    # Determine order and groups
    df = df_lsdd.copy()
    df[source_col] = df[source_col].astype(str).str.lower()
    sources = _sorted_sources(df[source_col].unique())

    data = [df.loc[df[source_col] == s, "LSDD"].dropna().values for s in sources]
    bp = plt.boxplot(data, labels=sources, vert=True, showfliers=False, patch_artist=True)

    # Colour boxes
    for patch, s in zip(bp["boxes"], sources):
        patch.set_facecolor(SOURCE_COLOURS.get(s, "#CCCCCC"))
        patch.set_edgecolor("#000000")
        patch.set_linewidth(1.4)

    for element in ["whiskers", "caps", "medians"]:
        for line in bp[element]:
            line.set_color("#000000")
            line.set_linewidth(1.4)

    plt.ylabel("LSDD")
    if title:
        plt.title(title)
    plt.grid(True, linestyle="--", alpha=0.35, axis="y")
    plt.xticks(rotation=20, ha="right")
    plt.tight_layout()
    if outpath:
        plt.savefig(outpath, dpi=300)
    return fig
