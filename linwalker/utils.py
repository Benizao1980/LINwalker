from __future__ import annotations

from typing import Iterable, List, Optional


def infer_max_level_from_lincode(series) -> int:
    """Infer maximum LIN depth from a LINcode column.

    LIN codes are treated as prefixes separated by '_' (common PubMLST export format).
    We take the maximum number of components observed.
    """
    max_n = 0
    for v in series.dropna().astype(str):
        n = len(v.split('_'))
        if n > max_n:
            max_n = n
    return max_n or 17


def parse_thresholds(spec: Optional[str], max_level: int) -> List[int]:
    """Parse a threshold specification.

    Supported:
    - None -> 1..max_level
    - '1-17'
    - '1,2,3,5'
    - '5' -> [5]
    """
    if not spec:
        return list(range(1, max_level + 1))

    s = str(spec).strip()
    if '-' in s:
        a, b = s.split('-', 1)
        a = int(a.strip())
        b = int(b.strip())
        if a > b:
            a, b = b, a
        a = max(1, a)
        b = min(max_level, b)
        return list(range(a, b + 1))

    if ',' in s:
        out = []
        for part in s.split(','):
            part = part.strip()
            if not part:
                continue
            out.append(int(part))
        out = sorted({t for t in out if 1 <= t <= max_level})
        return out

    # single int
    t = int(s)
    if 1 <= t <= max_level:
        return [t]
    return []
