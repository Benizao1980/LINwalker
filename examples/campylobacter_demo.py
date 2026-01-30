"""LINwalker demo (Campylobacter)

Run:
  python examples/campylobacter_demo.py
"""

import pandas as pd
from linwalker.diversification import lin_diversification
from linwalker.plotting import plot_diversification

df = pd.read_csv("../PATHSAFE_LINwalker_min.tsv", sep="\t")

div = lin_diversification(
    df,
    lin_col="LINcode (C. jejuni / C. coli cgMLST v2)",
    group_col="source"
)

plot_diversification(div, title="LIN diversification by host reservoir")
