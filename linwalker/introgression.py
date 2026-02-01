from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional, Tuple

import numpy as np
import pandas as pd

from .utils import lin_prefix, parse_thresholds


@dataclass
class IntrogressionResult:
    mixed_species: pd.DataFrame
    lsdd: pd.DataFrame


def mixed_species_fraction(
    df: pd.DataFrame,
    *,
    lin_col: str = "lin_code",
    species_col: str = "species",
    thresholds: Optional[str] = None,
    max_level: int = 17,
) -> pd.DataFrame:
    """Compute fraction of LIN prefixes that contain >1 species at each LIN level."""
    ks = parse_thresholds(thresholds, max_level)

    out = []
    for k in ks:
        pref = df[lin_col].astype(str).map(lambda s: lin_prefix(s, k))
        tmp = df.assign(lin_prefix=pref)
        n_species_per_prefix = tmp.groupby("lin_prefix")[species_col].nunique(dropna=True)
        mixed = (n_species_per_prefix > 1).mean() if len(n_species_per_prefix) else 0.0
        out.append({"lin_level": k, "mixed_fraction": float(mixed)})

    return pd.DataFrame(out)


def lsdd_by_level(
    df: pd.DataFrame,
    *,
    lin_col: str = "lin_code",
    species_col: str = "species",
    thresholds: Optional[str] = None,
    max_level: int = 17,
) -> pd.DataFrame:
    """A simple lineage-specific divergence density (LSDD)-like summary.

    For each LIN level k, compute the distribution of species composition per LIN prefix.

    Outputs per level:
      - n_prefixes
      - frac_mixed
      - mean_species_per_prefix
    """
    ks = parse_thresholds(thresholds, max_level)
    out = []
    for k in ks:
        pref = df[lin_col].astype(str).map(lambda s: lin_prefix(s, k))
        tmp = df.assign(lin_prefix=pref)
        counts = tmp.groupby("lin_prefix")[species_col].nunique(dropna=True)
        if len(counts) == 0:
            out.append({"lin_level": k, "n_prefixes": 0, "frac_mixed": 0.0, "mean_species_per_prefix": 0.0})
            continue
        out.append(
            {
                "lin_level": k,
                "n_prefixes": int(len(counts)),
                "frac_mixed": float((counts > 1).mean()),
                "mean_species_per_prefix": float(counts.mean()),
            }
        )
    return pd.DataFrame(out)

