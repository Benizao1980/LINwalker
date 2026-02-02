from __future__ import annotations

import argparse
import logging
from pathlib import Path
from typing import List

import pandas as pd

from . import __version__
from .diversification import lin_diversification, lin_diversification_rarefied
from .introgression import mixed_species_fraction, lsdd_by_level
from .stcc import stcc_concordance
from .outbreak import outbreak_descriptives
from .tree_export import export_tree_metadata
from .prep import prep_pubmlst_export
from .utils import resolve_column
from .plotting import (
    plot_diversification_curve,
    plot_mixed_species,
    plot_lsdd,
    plot_stcc_concordance,
    plot_outbreak_top_clusters,
    plot_outbreak_epicurve,
)


def _ensure_outdirs(outdir: Path):
    """Create a stable output structure under *outdir*.

    We keep outputs in predictable locations so downstream workflows (and CI)
    can check expected artifacts.
    """
    outdir.mkdir(parents=True, exist_ok=True)

    plots = outdir / "plots"
    tables = outdir / "tables"
    logs = outdir / "logs"
    derived = outdir / "derived"

    for p in (plots, tables, logs, derived):
        p.mkdir(parents=True, exist_ok=True)

    return plots, tables, logs, derived


def _setup_logging(logfile: Path, *, verbose: bool = False) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s | %(levelname)s | %(message)s",
        handlers=[logging.FileHandler(logfile), logging.StreamHandler()],
    )


def _read_tsv(path: str) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", low_memory=False)


def cmd_prep(args: argparse.Namespace) -> None:
    outdir = Path(args.outdir)
    plots, tables, logs, derived = _ensure_outdirs(outdir)
    _setup_logging(logs / "prep.log", verbose=getattr(args, "verbose", False))
    logging.info(f"LINwalker v{__version__} | prep")

    # Prep outputs are analysis-ready derived tables.
    # Keep them in outdir/derived so downstream commands can refer to them.
    _ = prep_pubmlst_export(
        Path(args.input),
        outdir=derived,
        prefix=args.prefix,
        bin_sources=not args.no_bin_sources,
        lin_col=args.lin_col,
        source_col=args.source_col,
        species_col=args.species_col,
        sample_col=args.sample_col,
        st_col=args.st_col,
        cc_col=args.cc_col,
        country_col=args.country_col,
        date_col=args.date_col,
        keep_cgmlst_matrix=args.keep_cgmlst_matrix,
    )

    logging.info(f"Wrote prep outputs to: {derived}")


def cmd_diversify(args: argparse.Namespace) -> None:
    outdir = Path(args.outdir)
    plots, tables, logs, _derived = _ensure_outdirs(outdir)
    _setup_logging(logs / "diversify.log", verbose=args.verbose)

    df = pd.read_csv(args.input, sep="\t", low_memory=False)

    args.lin_col = resolve_column(
        df,
        args.lin_col,
        ["lin_code", "LINcode", "LIN_code", "LIN", "lincode", "lin"],
    )
    if args.lin_col is None:
        raise ValueError(
            "Could not find a LIN column. Use --lin-col or ensure one of these exists: "
            "lin_code, LINcode, LIN_code, LIN, lincode, lin"
        )
    args.lin_col = resolve_column(
        df,
        args.lin_col,
        ["lin_code", "LINcode", "LIN_code", "LIN", "lincode", "lin"],
    )
    if args.lin_col is None:
        raise ValueError(
            "Missing LIN column. Tried --lin-col and common candidates (lin_code, LINcode, LIN_code, LIN, lincode, lin)."
        )
    args.lin_col = resolve_column(
        df,
        args.lin_col,
        ["lin_code", "LINcode", "LIN_code", "LIN", "lincode", "lin"],
    )
    if args.lin_col is None:
        raise ValueError(
            "Could not find a LIN code column. Use --lin-col to specify it. "
            f"Columns detected: {list(df.columns)}"
        )

    # Handle common PubMLST column naming variants (and older LINwalker outputs)
    args.lin_col = resolve_column(df, args.lin_col, [
        "lin_code", "LINcode", "LIN_code", "LIN", "lin", "lincode", "LINCODE",
    ])
    args.group_col = resolve_column(df, args.group_col, [
        "source", "Source", "group", "Group", "reservoir", "host", "metadata_source",
    ])
    if args.lin_col is None:
        raise ValueError(
            "Could not find a LIN column. Try --lin-col and check your input columns."
        )
    if args.group_col is None:
        raise ValueError(
            "Could not find a grouping/source column. Try --group-col and check your input columns."
        )
    div = lin_diversification(
        df,
        lin_col=args.lin_col,
        group_col=args.group_col,
        thresholds=args.thresholds,
        max_level=args.max_level,
    )

    div.table.to_csv(tables / "diversification.tsv", sep="\t", index=False)

    # Option A (default): reservoir-only plot
    reservoir = ["chicken", "ruminant", "pig", "wild bird"]
    df_res = div.table[div.table["source"].isin(reservoir)].copy()
    if not df_res.empty:
        plot_diversification_curve(
            df_res,
            outdir=plots,
            filename="diversification",
            formats=args.formats,
            title=args.title,
        )

    # Option B: all sources present
    plot_diversification_curve(
        div.table,
        outdir=plots,
        filename="diversification_all_sources",
        formats=args.formats,
        title=(args.title + " (all sources)") if args.title else "Diversification (all sources)",
    )

    # Rarefaction: normalise sample size per source to reduce sampling bias.
    if not args.no_rarefy:
        def _min_n(tbl, sources):
            sub = tbl[tbl["source"].isin(sources)]["source"].value_counts()
            if sub.empty:
                return None
            return int(sub.min())

        # Reservoir-only (default ecological view)
        if args.rarefy_n is not None:
            n_res = args.rarefy_n
        else:
            n_res = _min_n(df, reservoir)
        if n_res is not None and n_res >= 2:
            rare_res = lin_diversification_rarefied(
                df[df[args.group_col].isin(reservoir)].copy(),
                lin_col=args.lin_col,
                group_col=args.group_col,
                thresholds=args.thresholds,
                max_level=args.max_level,
                n_per_group=n_res,
                n_reps=args.rarefy_reps,
                seed=args.rarefy_seed,
            )
            rare_res.table.to_csv(tables / "diversification_rarefied.tsv", sep="\t", index=False)
            plot_diversification_curve(
                rare_res.table,
                outdir=plots,
                filename="diversification_rarefied",
                formats=args.formats,
                title=(
                    (args.title + f" (rarefied; n={n_res}, reps={args.rarefy_reps})")
                    if args.title
                    else f"Diversification (rarefied; n={n_res}, reps={args.rarefy_reps})"
                ),
            )

        # All sources (epidemiological view)
        if args.rarefy_n is not None:
            n_all = args.rarefy_n
        else:
            n_all = _min_n(df, df[args.group_col].dropna().unique())
        if n_all is not None and n_all >= 2:
            rare_all = lin_diversification_rarefied(
                df,
                lin_col=args.lin_col,
                group_col=args.group_col,
                thresholds=args.thresholds,
                max_level=args.max_level,
                n_per_group=n_all,
                n_reps=args.rarefy_reps,
                seed=args.rarefy_seed,
            )
            rare_all.table.to_csv(tables / "diversification_all_sources_rarefied.tsv", sep="\t", index=False)
            plot_diversification_curve(
                rare_all.table,
                outdir=plots,
                filename="diversification_all_sources_rarefied",
                formats=args.formats,
                title=(
                    (args.title + f" (all sources; rarefied; n={n_all}, reps={args.rarefy_reps})")
                    if args.title
                    else f"Diversification (all sources; rarefied; n={n_all}, reps={args.rarefy_reps})"
                ),
            )

    logging.info("[LINwalker] Wrote diversification outputs to: %s", outdir)

def cmd_introgress(args: argparse.Namespace) -> None:
    outdir = Path(args.outdir)
    plots, tables, logs, _derived = _ensure_outdirs(outdir)
    _setup_logging(logs / "introgress.log", verbose=args.verbose)

    df = pd.read_csv(args.input, sep="\t", low_memory=False)

    args.lin_col = resolve_column(
        df,
        args.lin_col,
        ["lin_code", "LINcode", "LIN_code", "LIN", "lincode", "lin"],
    )
    args.species_col = resolve_column(
        df,
        args.species_col,
        ["species", "Species", "organism", "Organism"],
    )
    if args.lin_col is None:
        raise ValueError("Could not find a LIN code column. Use --lin-col to specify it.")
    if args.species_col is None:
        raise ValueError("Could not find a species column. Use --species-col to specify it.")

    mix = mixed_species_fraction(
        df,
        lin_col=args.lin_col,
        species_col=args.species_col,
        thresholds=args.thresholds,
        max_level=args.max_level,
    )
    lsdd = lsdd_by_level(
        df,
        lin_col=args.lin_col,
        species_col=args.species_col,
        thresholds=args.thresholds,
        max_level=args.max_level,
    )

    # plotting.plot_lsdd expects an "lsdd" column
    if "lsdd" not in lsdd.columns and "mean_species_per_prefix" in lsdd.columns:
        lsdd = lsdd.rename(columns={"mean_species_per_prefix": "lsdd"})

    mix.to_csv(tables / "mixed_species.tsv", sep="\t", index=False)
    lsdd.to_csv(tables / "lsdd.tsv", sep="\t", index=False)

    plot_mixed_species(mix, outdir=plots, filename="mixed_species", formats=args.formats, title=args.title_mixed, max_level=args.max_level)
    plot_lsdd(lsdd, outdir=plots, filename="lsdd", formats=args.formats, title=args.title_lsdd)

    logging.info("Wrote introgression outputs to: %s", outdir)

def cmd_stcc(args: argparse.Namespace) -> None:
    outdir = Path(args.outdir)
    plots, tables, logs, _derived = _ensure_outdirs(outdir)
    _setup_logging(logs / "stcc.log", verbose=args.verbose)

    df = pd.read_csv(args.input, sep="\t", low_memory=False)

    args.lin_col = resolve_column(
        df,
        args.lin_col,
        ["lin_code", "LINcode", "LIN_code", "LIN", "lincode", "lin"],
    )
    if args.lin_col is None:
        raise ValueError(
            "Could not find a LIN column. Tried: lin_code, LINcode, LIN_code, LIN, lincode, lin."
        )

    res = stcc_concordance(
        df,
        lin_col=args.lin_col,
        st_col=args.st_col,
        cc_col=args.cc_col,
        thresholds=args.thresholds,
        max_level=args.max_level,
    )

    res.table.to_csv(tables / "stcc_concordance.tsv", sep="\t", index=False)
    plot_stcc_concordance(
        res.table,
        outdir=plots,
        filename="stcc_concordance",
        formats=args.formats,
        title=args.title,
    )
    logging.info("Wrote outputs to %s", outdir)

def cmd_tree(args: argparse.Namespace) -> None:
    outdir = Path(args.outdir)
    plots, tables, logs, _derived = _ensure_outdirs(outdir)
    _setup_logging(logs / "tree.log", verbose=getattr(args, "verbose", False))
    logging.info(f"LINwalker v{__version__} | tree")

    df = _read_tsv(args.input)

    args.lin_col = resolve_column(
        df,
        args.lin_col,
        ["lin_code", "LINcode", "LIN_code", "LIN", "lincode", "lin"],
    )
    if args.lin_col is None:
        raise ValueError("Could not find a LIN column. Use --lin-col to specify.")

    args.source_col = resolve_column(
        df,
        args.source_col,
        ["source", "reservoir", "host", "group"],
    )
    export_tree_metadata(
        df,
        outdir=tables,
        lin_col=args.lin_col,
        sample_col=args.sample_col,
        source_col=args.source_col,
        species_col=args.species_col,
        threshold=args.threshold,
        extra_cols=args.extra_cols,
    )
    logging.info(f"Wrote tree metadata to: {tables}")


def cmd_outbreak(args: argparse.Namespace) -> None:
    outdir = Path(args.outdir)
    plots, tables, logs, _derived = _ensure_outdirs(outdir)
    _setup_logging(logs / "outbreak.log", verbose=getattr(args, "verbose", False))
    logging.info(f"LINwalker v{__version__} | outbreak")

    df = _read_tsv(args.input)

    args.lin_col = resolve_column(
        df,
        args.lin_col,
        ["lin_code", "LINcode", "LIN_code", "LIN", "lincode", "lin"],
    )
    if args.lin_col is None:
        raise ValueError("Could not find a LIN column; try --lin-col.")

    args.source_col = resolve_column(df, args.source_col, ["source", "Source", "group", "reservoir", "host"])
    args.country_col = resolve_column(df, args.country_col, ["country", "Country", "location", "Location"])
    args.date_col = resolve_column(df, args.date_col, ["date", "collection_date", "isolation_date", "year", "Year"])

    res = outbreak_descriptives(
        df,
        lin_col=args.lin_col,
        source_col=args.source_col,
        country_col=args.country_col,
        date_col=args.date_col,
        top_n=args.top_n,
        top_threshold=args.top_threshold,
    )

    res.top_clusters.to_csv(tables / "top_clusters.tsv", sep="\t", index=False)
    res.cluster_source_counts.to_csv(tables / "cluster_source_counts.tsv", sep="\t", index=False)
    if res.cluster_country_counts is not None:
        res.cluster_country_counts.to_csv(tables / "cluster_country_counts.tsv", sep="\t", index=False)
    if res.cluster_date_counts is not None:
        res.cluster_date_counts.to_csv(tables / "cluster_date_counts.tsv", sep="\t", index=False)

    plot_outbreak_top_clusters(
        res.top_clusters,
        outdir=plots,
        filename="top_clusters",
        formats=args.formats,
        title=args.title,
    )
    if res.cluster_date_counts is not None:
        plot_outbreak_epicurve(
            res.cluster_date_counts,
            outdir=plots,
            filename="epicurve",
            formats=args.formats,
            title=(args.title + " (epi curve)" if args.title else None),
        )

    logging.info(f"Wrote outbreak outputs to: {outdir}")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="linwalker",
        description="LINwalker",
        allow_abbrev=False,
    )
    parser.add_argument("--verbose", action="store_true", help="Verbose logging")
    parser.add_argument("--version", action="version", version=f"LINwalker {__version__}")

    sub = parser.add_subparsers(dest="cmd", required=True)

    # prep
    p = sub.add_parser("prep", help="Prepare PubMLST export into analysis-ready tables", allow_abbrev=False)
    p.add_argument("--input", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--prefix", default="RUN", help="Prefix for output file basenames")
    p.add_argument("--no-bin-sources", action="store_true", help="Do not collapse/bin source labels")
    p.add_argument("--lin-col", default="LINcode")
    p.add_argument("--sample-col", default="isolate")
    p.add_argument("--source-col", default="source")
    p.add_argument("--species-col", default="species")
    p.add_argument("--st-col", default="ST")
    p.add_argument("--cc-col", default="clonal_complex")
    p.add_argument("--country-col", default="country")
    p.add_argument("--date-col", default="date")
    p.set_defaults(func=cmd_prep)

    # diversify
    p = sub.add_parser("diversify", help="Compute diversification curves and plot", allow_abbrev=False)
    p.add_argument("--input", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--lin-col", default="lin_code")
    p.add_argument("--group-col", default="source")
    p.add_argument("--thresholds", default=None, help='e.g. "1-17" or "1,3,5"')
    p.add_argument("--max-level", type=int, default=17)
    p.add_argument("--formats", nargs="+", default=["png", "svg"], help="Output image formats")
    p.add_argument("--title", default=None)
    # Rarefaction (sample-size normalisation)
    p.add_argument("--no-rarefy", action="store_true", help="Disable rarefaction curves")
    p.add_argument("--rarefy-reps", type=int, default=100, help="Rarefaction replicates")
    p.add_argument("--rarefy-seed", type=int, default=13, help="Random seed for rarefaction")
    p.add_argument(
        "--rarefy-n",
        type=int,
        default=None,
        help="Rarefaction depth (n per group). Default=min group size in plotted set",
    )
    p.set_defaults(func=cmd_diversify)

    # introgress
    p = sub.add_parser("introgress", help="Compute mixed-species curves and LSDD + plots", allow_abbrev=False)
    p.add_argument("--input", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--lin-col", default="lin_code")
    p.add_argument("--species-col", default="species")
    p.add_argument("--thresholds", default=None)
    p.add_argument("--max-level", type=int, default=17)
    p.add_argument("--formats", nargs="+", default=["png", "svg"])
    p.add_argument("--title-mixed", default=None)
    p.add_argument("--title-lsdd", default=None)
    p.set_defaults(func=cmd_introgress)

    # stcc
    p = sub.add_parser("stcc", help="Relate LIN thresholds to MLST ST and clonal complex", allow_abbrev=False)
    p.add_argument("--input", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--lin-col", default="lin_code")
    p.add_argument("--st-col", default="ST")
    p.add_argument("--cc-col", default="clonal_complex")
    p.add_argument("--thresholds", default=None)
    p.add_argument("--max-level", type=int, default=17)
    p.add_argument("--formats", nargs="+", default=["png", "svg"])
    p.add_argument("--title", default=None)
    p.set_defaults(func=cmd_stcc)

    # tree
    p = sub.add_parser("tree", help="Export Microreact/iTOL metadata to colour an existing tree", allow_abbrev=False)
    p.add_argument("--input", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--lin-col", default="lin_code")
    p.add_argument("--sample-col", default="isolate")
    p.add_argument("--source-col", default="source")
    p.add_argument("--species-col", default="species")
    p.add_argument("--threshold", type=int, default=12)
    p.add_argument("--extra-cols", default=None, help="Comma-separated extra metadata columns to include")
    p.set_defaults(func=cmd_tree)

    # outbreak
    p = sub.add_parser("outbreak", help="Descriptive outbreak/public-health summaries from LIN codes", allow_abbrev=False)
    p.add_argument("--input", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--lin-col", default="lin_code")
    p.add_argument("--source-col", default="source")
    p.add_argument("--species-col", default="species")
    p.add_argument("--country-col", default="country")
    p.add_argument("--date-col", default="date")
    p.add_argument("--top-n", type=int, default=25)
    p.add_argument("--top-threshold", type=int, default=12)
    p.add_argument(
        "--thresholds",
        default="1-17",
        help="Threshold range/string for plotting (accepted for compatibility; default 1-17)",
    )
    p.add_argument(
        "--max-level",
        type=int,
        default=17,
        help="Maximum LIN level to consider for plotting",
    )
    p.add_argument("--formats", nargs="+", default=["png", "svg"])
    p.add_argument("--title", default=None)
    p.set_defaults(func=cmd_outbreak)

    return parser


def main(argv: List[str] | None = None) -> None:
    parser = build_parser()
    args = parser.parse_args(argv)
    args.func(args)
