"""Command line interface for LINwalker."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from .prep import prepare_pubmlst_export
from .diversification import lin_diversification
from .introgression import mixed_species_summary, lsdd
from .stcc import stcc_concordance_by_threshold
from .plotting import plot_diversification_curve, plot_mixed_species, plot_lsdd_by_source, plot_stcc_concordance
from .tree_export import export_microreact_tsv, export_itol_colourstrip
from .utils import infer_max_level_from_lincode, parse_thresholds
from .outbreak import (
    per_isolate_cluster_sizes,
    summarise_levels,
    top_clusters_table,
    epi_curve,
    plot_cluster_size_summary,
    plot_cluster_size_boxplot,
)


def cmd_prep(args):
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    lin_df, meta_df, cg_df = prepare_pubmlst_export(
        args.input,
        outdir=outdir,
        prefix=args.prefix,
        source_col=args.source_col,
        species_col=args.species_col,
        bin_sources=not args.no_bin_sources,
    )

    # already written by prepare_pubmlst_export; return is mainly for API use
    print(f"[LINwalker] Wrote derived tables to: {outdir}")


def cmd_diversify(args):
    df = pd.read_csv(args.input, sep="\t", low_memory=False)
    inferred = infer_max_level_from_lincode(df[args.lin_col])
    max_level = min(int(args.max_level), int(inferred))
    thresholds = parse_thresholds(args.thresholds, max_level=max_level)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    div = lin_diversification(df, lin_col=args.lin_col, group_col=args.group_col, thresholds=thresholds, max_level=max_level)
    div.to_csv(outdir / "diversification.tsv", sep="\t", index=False)

    formats = args.formats

    # Reservoir-only plot (default ecological view)
    reservoir_sources = ["chicken", "ruminant", "pig", "wild bird"]
    div_res = div[div["source"].isin(reservoir_sources)]
    plot_diversification_curve(
        div_res,
        outpath_base=outdir / "diversification",
        title=args.title_reservoir or "Diversification by LIN threshold (reservoirs)",
        max_level=max_level,
        formats=formats,
    )

    # All sources (epi view)
    plot_diversification_curve(
        div,
        outpath_base=outdir / "diversification_all_sources",
        title=args.title_all or "Diversification by LIN threshold (all sources)",
        max_level=max_level,
        formats=formats,
    )

    print(f"[LINwalker] Wrote diversification outputs to: {outdir}")


def cmd_introgress(args):
    df = pd.read_csv(args.input, sep="\t", low_memory=False)
    inferred = infer_max_level_from_lincode(df[args.lin_col])
    max_level = min(int(args.max_level), int(inferred))
    thresholds = parse_thresholds(args.thresholds, max_level=max_level)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    mix = mixed_species_summary(df, lin_col=args.lin_col, species_col=args.species_col, thresholds=thresholds)
    mix.to_csv(outdir / "mixed_species.tsv", sep="\t", index=False)

    formats = args.formats
    plot_mixed_species(
        mix,
        outpath_base=outdir / "mixed_species",
        title=args.title_mixed or "Mixed-species LIN clusters vs LIN threshold",
        max_level=max_level,
        formats=formats,
    )

    lin_df_lsdd = lsdd(df, lin_col=args.lin_col, species_col=args.species_col, thresholds=thresholds)
    lin_df_lsdd.to_csv(outdir / "lsdd.tsv", sep="\t", index=False)

    plot_lsdd_by_source(
        lin_df_lsdd,
        outpath_base=outdir / "lsdd_by_source",
        title=args.title_lsdd or "LSDD by source",
        formats=formats,
        show_points=not args.no_points,
        style=args.lsdd_style,
    )

    if args.export_bridges:
        # Early discordance isolates
        early = lin_df_lsdd.sort_values(["LSDD", args.source_col, args.species_col]).head(args.n_bridges)
        early.to_csv(outdir / "bridge_isolates_early_lsdd.tsv", sep="\t", index=False)

        # Mixed clusters at a chosen threshold (default: minimum threshold)
        k = args.bridge_threshold if args.bridge_threshold is not None else min(thresholds)
        prefix = df[args.lin_col].astype(str).str.split("_").str[:k].str.join("_")
        tmp = df.copy()
        tmp["LIN_prefix"] = prefix
        # clusters where >1 species
        g = tmp.groupby("LIN_prefix")[args.species_col].nunique().reset_index(name="n_species")
        mixed = tmp.merge(g[g["n_species"] > 1], on="LIN_prefix")
        mixed.to_csv(outdir / f"mixed_species_clusters_LIN{k}.tsv", sep="\t", index=False)

    print(f"[LINwalker] Wrote introgression outputs to: {outdir}")


def cmd_stcc(args):
    df = pd.read_csv(args.input, sep="\t", low_memory=False)
    inferred = infer_max_level_from_lincode(df[args.lin_col])
    max_level = min(int(args.max_level), int(inferred))
    thresholds = parse_thresholds(args.thresholds, max_level=max_level)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    stcc = stcc_concordance_by_threshold(
        df,
        lin_col=args.lin_col,
        st_col=args.st_col,
        cc_col=args.cc_col,
        thresholds=thresholds,
        min_cluster_size=args.min_cluster_size,
    )

    stcc.to_csv(outdir / "stcc_concordance.tsv", sep="\t", index=False)

    plot_stcc_concordance(
        stcc,
        outpath_base=outdir / "stcc_concordance",
        title=args.title or "LIN ↔ MLST concordance vs LIN threshold",
        max_level=max(thresholds),
        formats=args.formats,
    )

    print(f"[LINwalker] Wrote ST/CC concordance outputs to: {outdir}")


def cmd_tree(args):
    df = pd.read_csv(args.input, sep="\t", low_memory=False)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    micro = export_microreact_tsv(
        df,
        outdir / "microreact.tsv",
        id_col=args.id_col,
        source_col=args.source_col,
        species_col=args.species_col,
        lin_col=args.lin_col,
    )

    itol = export_itol_colourstrip(
        df,
        outdir / f"itol_colourstrip_{args.category}.txt",
        id_col=args.id_col,
        category_col=args.category,
        title=args.title or "LINwalker",
    )

    print(f"[LINwalker] Wrote tree metadata: {micro} | {itol}")


def cmd_outbreak(args):
    df = pd.read_csv(args.input, sep="\t", low_memory=False)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    inferred = infer_max_level_from_lincode(df[args.lin_col])
    max_level = min(int(args.max_level), int(inferred))
    thresholds = parse_thresholds(args.thresholds, max_level=max_level)

    # Per-isolate cluster sizes across levels
    per_iso = per_isolate_cluster_sizes(df, lin_col=args.lin_col, thresholds=thresholds)
    per_iso.to_csv(outdir / "per_isolate_cluster_sizes.tsv", sep="\t", index=False)

    # Level summaries
    levels = summarise_levels(df, lin_col=args.lin_col, thresholds=thresholds)
    levels.to_csv(outdir / "cluster_levels_summary.tsv", sep="\t", index=False)

    # Top clusters at a chosen threshold (default: max threshold)
    top_thr = int(args.top_threshold) if args.top_threshold is not None else int(max(thresholds))
    top_thr = min(top_thr, int(max(thresholds)))
    top = top_clusters_table(
        df,
        lin_col=args.lin_col,
        threshold=top_thr,
        n=int(args.top_n),
        source_col=args.source_col,
        species_col=args.species_col,
        country_col=args.country_col,
        id_col=args.id_col,
    )
    top.to_csv(outdir / f"top_clusters_LIN{top_thr}.tsv", sep="\t", index=False)

    # Optional epi curve
    if args.date_col:
        try:
            ec = epi_curve(df, date_col=args.date_col, lin_col=args.lin_col, threshold=top_thr)
            ec.to_csv(outdir / f"epi_curve_LIN{top_thr}.tsv", sep="\t", index=False)
        except Exception as e:
            print(f"[LINwalker] Warning: could not compute epi curve ({e})")

    # Plots
    plot_cluster_size_summary(
        levels,
        outpath_base=outdir / "cluster_size_summary",
        title=args.title or "Per-isolate LIN cluster size vs LIN threshold",
        formats=args.formats,
        max_level=max(thresholds),
    )
    plot_cluster_size_boxplot(
        per_iso,
        outpath_base=outdir / "cluster_size_boxplot",
        title=args.title_box or "Per-isolate cluster sizes (box + points)",
        formats=args.formats,
        max_level=max(thresholds),
        sample_points=int(args.sample_points),
    )

    print(f"[LINwalker] Wrote outbreak/public-health outputs to: {outdir}")


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="linwalker", description="LINwalker")
    sub = p.add_subparsers(dest="cmd", required=True)

    # prep
    s = sub.add_parser("prep", help="Prepare PubMLST export into analysis-ready tables")
    s.add_argument("--input", required=True, help="PubMLST export TSV/TSV.GZ")
    s.add_argument("--outdir", required=True, help="Output directory for derived tables")
    s.add_argument("--prefix", default=None, help="Output filename prefix")
    s.add_argument("--source-col", default="source", help="Source column")
    s.add_argument("--species-col", default="species", help="Species column")
    s.add_argument("--no-bin-sources", action="store_true", help="Do not collapse source labels into stable bins")
    s.set_defaults(func=cmd_prep)

    # diversify
    s = sub.add_parser("diversify", help="Compute diversification curves and plot")
    s.add_argument("--input", required=True, help="LINwalker minimal TSV")
    s.add_argument("--lin-col", default="LINcode")
    s.add_argument("--group-col", default="source")
    s.add_argument("--thresholds", default=None, help="Thresholds, e.g. '1-17' or '1,2,5'")
    s.add_argument("--max-level", type=int, default=17, help="Maximum LIN level to analyze/plot (default: 17)")
    s.add_argument("--formats", nargs="+", default=["png", "svg"], help="Output formats (png svg)")
    s.add_argument("--outdir", required=True)
    s.add_argument("--title-reservoir", default=None)
    s.add_argument("--title-all", default=None)
    s.set_defaults(func=cmd_diversify)

    # introgress
    s = sub.add_parser("introgress", help="Compute mixed-species curves and LSDD + plots")
    s.add_argument("--input", required=True, help="LINwalker minimal TSV")
    s.add_argument("--lin-col", default="LINcode")
    s.add_argument("--species-col", default="species")
    s.add_argument("--source-col", default="source")
    s.add_argument("--thresholds", default=None, help="Thresholds, e.g. '1-17' or '1,2,5'")
    s.add_argument("--max-level", type=int, default=17, help="Maximum LIN level to analyze/plot (default: 17)")
    s.add_argument("--formats", nargs="+", default=["png", "svg"], help="Output formats (png svg)")
    s.add_argument("--outdir", required=True)
    s.add_argument("--title-mixed", default=None)
    s.add_argument("--title-lsdd", default=None)
    s.add_argument("--no-points", action="store_true", help="Do not overlay per-isolate points on LSDD plot")
    s.add_argument("--lsdd-style", default="violin_box", choices=["violin_box", "violin", "box"], help="LSDD plot style")
    s.add_argument("--export-bridges", action="store_true", help="Write bridge isolate / mixed-cluster tables")
    s.add_argument("--n-bridges", type=int, default=250, help="Number of early-LSDD isolates to export")
    s.add_argument("--bridge-threshold", type=int, default=None, help="LIN threshold for exporting mixed clusters")
    s.set_defaults(func=cmd_introgress)

    # stcc
    s = sub.add_parser("stcc", help="Relate LIN thresholds to MLST ST and clonal complex")
    s.add_argument("--input", required=True, help="Metadata TSV containing LINcode, ST and clonal complex")
    s.add_argument("--lin-col", default="LINcode")
    s.add_argument("--st-col", default="ST (MLST)")
    s.add_argument("--cc-col", default="clonal_complex (MLST)")
    s.add_argument("--thresholds", default=None, help="Thresholds, e.g. '1-17' or '1,2,5'")
    s.add_argument("--max-level", type=int, default=17, help="Maximum LIN level to analyze/plot (default: 17)")
    s.add_argument("--formats", nargs="+", default=["png", "svg"], help="Output formats (png svg)")
    s.add_argument("--min-cluster-size", type=int, default=1, help="Exclude LIN clusters smaller than this size")
    s.add_argument("--outdir", required=True)
    s.add_argument("--title", default=None)
    s.set_defaults(func=cmd_stcc)

    # tree exports
    s = sub.add_parser("tree", help="Export Microreact/iTOL metadata to colour an existing tree")
    s.add_argument("--input", required=True, help="LINwalker minimal TSV")
    s.add_argument("--outdir", required=True)
    s.add_argument("--id-col", default="isolate")
    s.add_argument("--lin-col", default="LINcode")
    s.add_argument("--source-col", default="source")
    s.add_argument("--species-col", default="species")
    s.add_argument("--category", default="source", help="Column to colour-by for iTOL (default: source)")
    s.add_argument("--title", default=None, help="Label used in iTOL dataset")
    s.set_defaults(func=cmd_tree)

    # outbreak / public health summaries
    s = sub.add_parser("outbreak", help="Descriptive outbreak/public-health summaries from LIN codes")
    s.add_argument("--input", required=True, help="LINwalker minimal TSV")
    s.add_argument("--lin-col", default="LINcode", help="Full LIN code column")
    s.add_argument("--source-col", default="source", help="Source column")
    s.add_argument("--species-col", default="species", help="Species column")
    s.add_argument("--country-col", default="country", help="Country/region column")
    s.add_argument("--date-col", default=None, help="Optional collection date column")
    s.add_argument("--top-n", type=int, default=10, help="Number of top clusters to export")
    s.add_argument("--top-threshold", type=int, default=12, help="LIN threshold for top-cluster table")
    s.add_argument("--max-level", type=int, default=17, help="Maximum LIN threshold to consider")
    s.add_argument("--formats", nargs="+", default=["png", "svg"], help="Output figure formats")
    s.add_argument("--title", default=None, help="Optional title prefix")
    s.add_argument("--outdir", required=True)
    s.set_defaults(func=cmd_outbreak)

    return p


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    args.func(args)


