# linwalker/cli.py

from __future__ import annotations
import argparse
from pathlib import Path
import pandas as pd

from .prep import prepare_pubmlst_export
from .diversification import lin_diversification
from .introgression import mixed_species_summary, lsdd
from .stcc import stcc_concordance_by_threshold
from .plotting import (
    plot_diversification,
    plot_mixed_species,
    plot_lsdd_by_source,
    plot_stcc_concordance,
)


def _ensure_outdir(outdir: str) -> Path:
    p = Path(outdir)
    p.mkdir(parents=True, exist_ok=True)
    return p


def cmd_prep(args):
    outdir = _ensure_outdir(args.outdir)
    prepare_pubmlst_export(
        args.input,
        outdir=str(outdir),
        prefix=args.prefix,
        source_col=args.source_col,
        species_col=args.species_col,
        bin_sources=not args.no_bin_sources,
    )
    print(f"[LINwalker] Wrote derived tables to: {outdir}")


def cmd_diversify(args):
    outdir = _ensure_outdir(args.outdir)
    df = pd.read_csv(args.input, sep="\t", low_memory=False)
    div = lin_diversification(df, lin_col=args.lin_col, group_col=args.group_col)
    (outdir / "diversification.tsv").write_text(div.to_csv(sep="\t", index=False))

    exclude = set() if args.include_human_other else None
    plot_diversification(div, title=args.title, outpath=str(outdir / "diversification.png"), exclude_sources=exclude)
    plot_diversification(div, title=args.title, outpath=str(outdir / "diversification.svg"), exclude_sources=exclude)

    # always export a full version as a convenience
    plot_diversification(div, title=f"{args.title} (all sources)", outpath=str(outdir / "diversification_all_sources.png"), exclude_sources=set())
    plot_diversification(div, title=f"{args.title} (all sources)", outpath=str(outdir / "diversification_all_sources.svg"), exclude_sources=set())

    print(f"[LINwalker] Wrote diversification outputs to: {outdir}")


def cmd_introgress(args):
    outdir = _ensure_outdir(args.outdir)
    df = pd.read_csv(args.input, sep="\t", low_memory=False)

    mix = mixed_species_summary(df, lin_col=args.lin_col, species_col=args.species_col)
    (outdir / "mixed_species.tsv").write_text(mix.to_csv(sep="\t", index=False))
    plot_mixed_species(mix, title=args.title_mixed, outpath=str(outdir / "mixed_species.png"))
    plot_mixed_species(mix, title=args.title_mixed, outpath=str(outdir / "mixed_species.svg"))

    ls = lsdd(df, lin_col=args.lin_col, species_col=args.species_col)
    (outdir / "lsdd.tsv").write_text(ls.to_csv(sep="\t", index=False))
    plot_lsdd_by_source(ls, source_col=args.source_col, title=args.title_lsdd, outpath=str(outdir / "lsdd_by_source.png"))
    plot_lsdd_by_source(ls, source_col=args.source_col, title=args.title_lsdd, outpath=str(outdir / "lsdd_by_source.svg"))

    print(f"[LINwalker] Wrote introgression outputs to: {outdir}")


def cmd_stcc(args):
    outdir = _ensure_outdir(args.outdir)
    df = pd.read_csv(args.input, sep="\t", low_memory=False)

    stcc = stcc_concordance_by_threshold(
        df,
        lin_col=args.lin_col,
        st_col=args.st_col,
        cc_col=args.cc_col,
        min_cluster_size=args.min_cluster_size,
    )
    (outdir / "stcc_concordance.tsv").write_text(stcc.to_csv(sep="\t", index=False))

    plot_stcc_concordance(stcc, title=args.title, outpath=str(outdir / "stcc_concordance.png"))
    plot_stcc_concordance(stcc, title=args.title, outpath=str(outdir / "stcc_concordance.svg"))

    print(f"[LINwalker] Wrote ST/CC concordance outputs to: {outdir}")


def build_parser():
    p = argparse.ArgumentParser(prog="linwalker", description="LINwalker")
    sub = p.add_subparsers(dest="cmd", required=True)

    sp = sub.add_parser("prep", help="Prepare PubMLST export into analysis-ready tables")
    sp.add_argument("--input", required=True, help="PubMLST export TSV/TSV.GZ")
    sp.add_argument("--outdir", required=True, help="Output directory for derived tables")
    sp.add_argument("--prefix", default=None, help="Output filename prefix (default: derived from input filename)")
    sp.add_argument("--source-col", default="source")
    sp.add_argument("--species-col", default="species")
    sp.add_argument("--no-bin-sources", action="store_true", help="Do not collapse source labels into stable bins")
    sp.set_defaults(func=cmd_prep)

    sd = sub.add_parser("diversify", help="Compute diversification curves and plot")
    sd.add_argument("--input", required=True, help="LINwalker minimal TSV")
    sd.add_argument("--lin-col", default="LINcode", help="Full LIN code column (default: LINcode)")
    sd.add_argument("--group-col", default="source", help="Grouping column (default: source)")
    sd.add_argument("--outdir", required=True)
    sd.add_argument("--title", default="Unique LIN IDs vs. LIN threshold by source (Campylobacter cgMLST v2)")
    sd.add_argument("--include-human-other", action="store_true", help="Include human and other in diversification plot (may dominate y-axis)")
    sd.set_defaults(func=cmd_diversify)

    si = sub.add_parser("introgress", help="Compute mixed-species curves and LSDD + plots")
    si.add_argument("--input", required=True, help="LINwalker minimal TSV")
    si.add_argument("--lin-col", default="LINcode", help="Full LIN code column (default: LINcode)")
    si.add_argument("--species-col", default="species", help="Species column (default: species)")
    si.add_argument("--source-col", default="source", help="Source column (default: source)")
    si.add_argument("--outdir", required=True)
    si.add_argument("--title-mixed", default="Mixed-species LIN clusters vs. LIN threshold (Campylobacter cgMLST v2)")
    si.add_argument("--title-lsdd", default="LIN Species Discordance Depth (LSDD) by source")
    si.set_defaults(func=cmd_introgress)

    ss = sub.add_parser("stcc", help="Relate LIN thresholds to MLST ST and clonal complex")
    ss.add_argument("--input", required=True, help="Metadata TSV containing LINcode, ST and clonal complex")
    ss.add_argument("--lin-col", default="LINcode")
    ss.add_argument("--st-col", default="ST (MLST)")
    ss.add_argument("--cc-col", default="clonal_complex (MLST)")
    ss.add_argument("--min-cluster-size", type=int, default=2, help="Exclude LIN clusters smaller than this size")
    ss.add_argument("--outdir", required=True)
    ss.add_argument("--title", default="Concordance between LIN thresholds and MLST ST / clonal complex")
    ss.set_defaults(func=cmd_stcc)

    return p


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
