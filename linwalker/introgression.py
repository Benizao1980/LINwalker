"""Introgression / species-boundary erosion summaries.

Outputs:
- mixed-species LIN clusters vs threshold
- LSDD (LIN Species Discordance Depth) per isolate
- optional 'bridge' tables highlighting mixed clusters and early-discordance isolates
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Tuple

import pandas as pd

from .plotting import plot_lsdd_by_source, plot_mixed_species
from .utils import infer_max_level_from_lincode, parse_thresholds


def _prefix_at_k(series: pd.Series, k: int) -> pd.Series:
    parts = series.astype(str).str.split("_", expand=True)
    return parts.iloc[:, :k].apply(lambda r: "_".join(r.values.astype(str)), axis=1)


def mixed_species_summary(
    df: pd.DataFrame,
    lin_col: str = "LINcode",
    species_col: str = "species",
    thresholds: Optional[str] = None,
) -> pd.DataFrame:
    max_level = infer_max_level_from_lincode(df[lin_col])
    ks = parse_thresholds(thresholds, max_level)

    rows = []
    for k in ks:
        prefix = _prefix_at_k(df[lin_col], k)
        tmp = df[[species_col]].copy()
        tmp["prefix"] = prefix

        # how many unique species per cluster
        sp_counts = tmp.groupby("prefix")[species_col].nunique()
        n_clusters = int(sp_counts.shape[0])
        n_mixed = int((sp_counts > 1).sum())
        prop_mixed = (n_mixed / n_clusters) if n_clusters else 0.0

        rows.append(
            {
                "threshold": k,
                "n_clusters": n_clusters,
                "n_mixed_species_clusters": n_mixed,
                "prop_mixed_species_clusters": prop_mixed,
            }
        )

    return pd.DataFrame(rows)


def lsdd(
    df: pd.DataFrame,
    lin_col: str = "LINcode",
    species_col: str = "species",
    thresholds: Optional[str] = None,
) -> pd.DataFrame:
    """Compute per-isolate LSDD.

    LSDD is the earliest LIN threshold at which the isolate's assigned species differs
    from the majority species of its LIN cluster.

    If no discordance is observed, the isolate gets LSDD=max_level.
    """
    max_level = infer_max_level_from_lincode(df[lin_col])
    ks = parse_thresholds(thresholds, max_level)

    out = df.copy()
    out["LSDD"] = max(ks) if ks else max_level
    out["LSDD_majority_species"] = pd.NA
    out["LSDD_cluster_size"] = pd.NA
    out["LSDD_prefix"] = pd.NA

    # Precompute prefixes for each k to avoid repeated splits
    prefix_by_k = {k: _prefix_at_k(df[lin_col], k) for k in ks}

    for k in ks:
        prefix = prefix_by_k[k]
        tmp = pd.DataFrame({"prefix": prefix, species_col: df[species_col].astype(str)})

        # majority species per prefix
        grp = tmp.groupby("prefix")[species_col]
        majority = grp.agg(lambda s: s.value_counts().index[0])
        size = grp.size()

        maj_map = majority.to_dict()
        size_map = size.to_dict()

        assigned = df[species_col].astype(str)
        maj_species = prefix.map(maj_map)
        discordant = assigned.ne(maj_species)

        # only set where not set yet
        not_set = out["LSDD"].astype(int) == (max(ks) if ks else max_level)
        to_set = discordant & not_set
        if to_set.any():
            out.loc[to_set, "LSDD"] = k
            out.loc[to_set, "LSDD_majority_species"] = maj_species[to_set].values
            out.loc[to_set, "LSDD_cluster_size"] = prefix[to_set].map(size_map).values
            out.loc[to_set, "LSDD_prefix"] = prefix[to_set].values

    return out


def bridge_clusters(
    df: pd.DataFrame,
    lin_col: str = "LINcode",
    species_col: str = "species",
    thresholds: Optional[str] = None,
    min_cluster_size: int = 5,
) -> pd.DataFrame:
    """List mixed-species LIN clusters that are potentially 'bridging' species.

    Returns one row per (threshold, prefix) for clusters with >1 species.
    """
    max_level = infer_max_level_from_lincode(df[lin_col])
    ks = parse_thresholds(thresholds, max_level)

    rows = []
    for k in ks:
        prefix = _prefix_at_k(df[lin_col], k)
        tmp = df[[species_col]].copy()
        tmp["prefix"] = prefix

        grp = tmp.groupby("prefix")[species_col]
        n = grp.size()
        nsp = grp.nunique()

        mixed_prefixes = nsp[(nsp > 1) & (n >= min_cluster_size)].index
        for pfx in mixed_prefixes:
            counts = grp.get_group(pfx).value_counts().to_dict()
            rows.append(
                {
                    "threshold": k,
                    "prefix": pfx,
                    "cluster_size": int(n[pfx]),
                    "n_species": int(nsp[pfx]),
                    "species_counts": ";".join([f"{s}:{c}" for s, c in sorted(counts.items(), key=lambda x: (-x[1], x[0]))]),
                }
            )

    return pd.DataFrame(rows)


def run_introgression(
    inpath: str,
    outdir: str,
    lin_col: str = "LINcode",
    species_col: str = "species",
    source_col: str = "source",
    thresholds: Optional[str] = None,
    formats: str = "png,svg",
    title_mixed: Optional[str] = None,
    title_lsdd: Optional[str] = None,
    show_points: bool = True,
    write_bridges: bool = True,
    min_bridge_size: int = 5,
) -> None:
    outdir_p = Path(outdir)
    outdir_p.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(inpath, sep="\t", low_memory=False)

    mix = mixed_species_summary(df, lin_col=lin_col, species_col=species_col, thresholds=thresholds)
    mix.to_csv(outdir_p / "mixed_species_summary.tsv", sep="\t", index=False)

    max_level = int(mix["threshold"].max()) if not mix.empty else infer_max_level_from_lincode(df[lin_col])
    fmts = [f.strip() for f in formats.split(',') if f.strip()]

    plot_mixed_species(
        mix,
        max_level=max_level,
        title=title_mixed or "Mixed-species LIN clusters by threshold",
        outpath=outdir_p / "mixed_species",
        formats=fmts,
    )

    lsdd_df = lsdd(df, lin_col=lin_col, species_col=species_col, thresholds=thresholds)
    lsdd_df.to_csv(outdir_p / "lsdd_per_isolate.tsv", sep="\t", index=False)

    plot_lsdd_by_source(
        lsdd_df,
        source_col=source_col,
        max_level=max_level,
        title=title_lsdd or "LSDD by source",
        outpath=outdir_p / "lsdd_by_source",
        formats=fmts,
        show_points=show_points,
    )

    if write_bridges:
        bridges = bridge_clusters(
            df,
            lin_col=lin_col,
            species_col=species_col,
            thresholds=thresholds,
            min_cluster_size=min_bridge_size,
        )
        bridges.to_csv(outdir_p / "bridge_mixed_species_clusters.tsv", sep="\t", index=False)

        # a quick 'early discordance' isolate list
        if "id" in lsdd_df.columns:
            early = lsdd_df.sort_values("LSDD").head(200)
            early.to_csv(outdir_p / "bridge_isolates_early_lsdd.tsv", sep="\t", index=False)

