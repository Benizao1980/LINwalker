from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional

import pandas as pd

from .utils import lin_prefix, parse_thresholds


@dataclass
class STCCResult:
    table: pd.DataFrame


def stcc_concordance(
    df: pd.DataFrame,
    *,
    lin_col: str = "lin_code",
    st_col: str = "ST",
    cc_col: str = "clonal_complex",
    thresholds: Optional[str] = None,
    max_level: int = 17,
) -> STCCResult:
    """At each LIN level, how well do LIN clusters correspond to ST / CC?

    Outputs per lin_level:
      - n_isolates
      - n_lin_clusters
      - n_ST
      - n_CC
      - mean_purity_ST (mean of max-ST fraction per LIN cluster)
      - mean_purity_CC
    """
    ks = parse_thresholds(thresholds, max_level)

    # Keep only rows with LIN and at least one of ST/CC
    x = df.copy()
    x = x[x[lin_col].notna()]

    records = []
    for k in ks:
        prefix = x[lin_col].map(lambda s: lin_prefix(s, k))
        tmp = x.assign(lin_prefix=prefix)

        # Purity relative to ST
        st_purity = None
        if st_col in tmp.columns:
            st_counts = tmp.groupby(["lin_prefix", st_col], dropna=False).size().reset_index(name="n")
            st_max = st_counts.groupby("lin_prefix")["n"].max()
            st_total = tmp.groupby("lin_prefix").size()
            st_purity = (st_max / st_total).mean()

        # Purity relative to CC
        cc_purity = None
        if cc_col in tmp.columns:
            cc_counts = tmp.groupby(["lin_prefix", cc_col], dropna=False).size().reset_index(name="n")
            cc_max = cc_counts.groupby("lin_prefix")["n"].max()
            cc_total = tmp.groupby("lin_prefix").size()
            cc_purity = (cc_max / cc_total).mean()

        rec = {
            "lin_level": k,
            "n_isolates": int(len(tmp)),
            "n_lin_clusters": int(tmp["lin_prefix"].nunique()),
            "n_ST": int(tmp[st_col].nunique()) if st_col in tmp.columns else pd.NA,
            "n_CC": int(tmp[cc_col].nunique()) if cc_col in tmp.columns else pd.NA,
            "mean_purity_ST": float(st_purity) if st_purity is not None else pd.NA,
            "mean_purity_CC": float(cc_purity) if cc_purity is not None else pd.NA,
        }
        records.append(rec)

    out = pd.DataFrame.from_records(records).sort_values("lin_level")
    return STCCResult(table=out)
