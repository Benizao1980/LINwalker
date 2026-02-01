from __future__ import annotations

import re
from typing import Iterable, List, Optional


def parse_thresholds(thresholds: Optional[str], max_level: int) -> List[int]:
    """Parse a thresholds spec into an ordered list of LIN levels.

    Examples
    --------
    "1-17" -> [1,2,...,17]
    "1,3,5" -> [1,3,5]
    None -> [1..max_level]

    Notes
    -----
    LIN levels are 1-indexed.
    """
    if thresholds is None or str(thresholds).strip() == "":
        return list(range(1, max_level + 1))

    s = str(thresholds).strip()
    out: List[int] = []

    # split on commas/spaces
    for part in re.split(r"[\s,]+", s):
        if not part:
            continue
        if "-" in part:
            a, b = part.split("-", 1)
            a_i = int(a)
            b_i = int(b)
            if a_i <= b_i:
                out.extend(list(range(a_i, b_i + 1)))
            else:
                out.extend(list(range(a_i, b_i - 1, -1)))
        else:
            out.append(int(part))

    # normalize
    out = [k for k in out if 1 <= k <= max_level]
    # unique but keep order
    seen = set()
    uniq: List[int] = []
    for k in out:
        if k not in seen:
            uniq.append(k)
            seen.add(k)
    return uniq


def lin_prefix(lin: str, level: int) -> str:
    """Return LIN prefix at a given 1-indexed level.

    Accepts LIN codes like "0_0_0_1" or "0.0.0.1".
    """
    if lin is None:
        return ""
    s = str(lin).strip()
    if s == "":
        return ""
    if "_" in s:
        parts = s.split("_")
    elif "." in s:
        parts = s.split(".")
    else:
        # fallback: split on whitespace
        parts = re.split(r"\s+", s)

    # remove empty tokens
    parts = [p for p in parts if p != ""]
    if level <= 0:
        return ""
    return "_".join(parts[:level])


def coerce_source(val: str) -> str:
    if val is None:
        return "other"
    s = str(val).strip().lower()
    if s in {"human", "clinical", "patient"}:
        return "human"
    if s in {"chicken", "poultry"}:
        return "chicken"
    if s in {"cattle", "cow", "sheep", "goat", "ruminant"}:
        return "ruminant"
    if s in {"pig", "swine"}:
        return "pig"
    if s in {"wildbird", "wild bird", "bird", "avian"}:
        return "wild bird"
    return s if s else "other"


def resolve_column(df, preferred: str | None, candidates: list[str]) -> str | None:
    """Resolve a column name from a dataframe.

    If *preferred* exists in the dataframe, it is returned.
    Otherwise try *candidates* in order.
    Returns None if no column can be resolved.
    """
    if preferred and preferred in df.columns:
        return preferred
    for c in candidates:
        if c in df.columns:
            return c
    return None
