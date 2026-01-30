"""LINwalker demo (Campylobacter / PubMLST)

Run:
  python -m linwalker --help
  python -m linwalker prep --input PATHSAFE_pubmlst_export.tsv.gz --outdir data/derived --prefix PATHSAFE
  python -m linwalker diversify --input data/derived/PATHSAFE_LINwalker_min.tsv --lin-col LINcode --group-col source --outdir results/diversification
  python -m linwalker introgress --input data/derived/PATHSAFE_LINwalker_min.tsv --lin-col LINcode --species-col species --outdir results/introgression
"""
