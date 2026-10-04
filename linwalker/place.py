from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from pathlib import Path
import re
from typing import Iterable, Sequence

import numpy as np
import pandas as pd

from .utils import lin_prefix


CAMPY_V2_LIN_DIFF_THRESHOLDS = (
    1119, 1085, 982, 914, 857, 680, 445, 343, 183,
    86, 43, 10, 7, 5, 3, 2, 1, 0,
)

CAMPY_V2_CGC_THRESHOLDS = (200, 100, 50, 25, 10, 5)

DEFAULT_MISSING_TOKENS = {
    "", "nan", "none", "null", "na", "n/a", "<na>", "-", "x", "i",
}


@dataclass
class PlacementResult:
    summary: pd.DataFrame
    by_threshold: pd.DataFrame
    cgc2: pd.DataFrame
    pairwise: pd.DataFrame
    label_map: pd.DataFrame


def parse_diff_thresholds(value: str | Sequence[int] | None) -> list[int]:
    if value is None:
        return list(CAMPY_V2_LIN_DIFF_THRESHOLDS)
    if isinstance(value, str):
        parts = [x for x in re.split(r"[;,\s]+", value.strip()) if x]
        out = [int(x) for x in parts]
    else:
        out = [int(x) for x in value]
    if not out:
        raise ValueError("At least one LIN difference threshold is required.")
    if any(x < 0 for x in out):
        raise ValueError("LIN difference thresholds must be >=0.")
    if any(a < b for a, b in zip(out, out[1:])):
        raise ValueError("LIN difference thresholds must be ordered coarse-to-fine (descending).")
    return out


def _normalise_allele(value, missing_tokens: set[str]):
    if pd.isna(value):
        return None
    s = str(value).strip()
    if s.lower() in missing_tokens:
        return None
    # Excel commonly turns integer allele IDs into strings such as '12.0'.
    if re.fullmatch(r"\d+\.0", s):
        s = s[:-2]
    return s


def read_profiles(
    path: str | Path,
    *,
    sheet: str = "all",
    sample_id_col: str | None = None,
    locus_prefix: str = "CAMP",
) -> pd.DataFrame:
    """Read sample x locus allele profiles.

    For a PubMLST Genome Comparator XLSX file, the 'all' worksheet has loci in
    rows and genomes in columns; this function transposes it automatically.

    TSV/CSV input is expected to have genomes in rows. Locus columns are
    identified by *locus_prefix*.
    """
    path = Path(path)
    suffix = path.suffix.lower()

    if suffix in {".xlsx", ".xlsm", ".xls"}:
        raw = pd.read_excel(path, sheet_name=sheet, dtype=str)
        if raw.empty:
            raise ValueError(f"Worksheet '{sheet}' is empty: {path}")
        locus_col = "Locus" if "Locus" in raw.columns else raw.columns[0]
        known_meta = {
            locus_col,
            "Full name",
            "Product",
            "Sequence length",
            " Genome position",
            "Genome position",
            "Reference genome",
        }
        sample_cols = [c for c in raw.columns if c not in known_meta]
        raw = raw[raw[locus_col].notna()].copy()
        raw[locus_col] = raw[locus_col].astype(str)
        raw = raw[raw[locus_col].str.startswith(locus_prefix)]
        if raw.empty:
            raise ValueError(
                f"No loci beginning '{locus_prefix}' found in Genome Comparator sheet '{sheet}'."
            )
        if raw[locus_col].duplicated().any():
            dup = raw.loc[raw[locus_col].duplicated(), locus_col].iloc[0]
            raise ValueError(f"Duplicate locus in profile worksheet: {dup}")
        prof = raw.set_index(locus_col)[sample_cols].T
        prof.index = prof.index.map(str)
        prof.index.name = "sample_id"
        return prof

    sep = "," if suffix == ".csv" else "\t"
    raw = pd.read_csv(path, sep=sep, dtype=str, low_memory=False)
    if raw.empty:
        raise ValueError(f"Profile table is empty: {path}")
    if sample_id_col is None:
        sample_id_col = raw.columns[0]
    if sample_id_col not in raw.columns:
        raise ValueError(f"Sample ID column not found: {sample_id_col}")
    loci = [c for c in raw.columns if str(c).startswith(locus_prefix)]
    if not loci:
        raise ValueError(f"No locus columns beginning '{locus_prefix}' found in {path}")
    prof = raw.set_index(sample_id_col)[loci].copy()
    prof.index = prof.index.map(str)
    prof.index.name = "sample_id"
    return prof


def _label_tokens(label: str) -> set[str]:
    s = str(label).strip()
    tokens = {s}
    stem = Path(s).stem
    tokens.add(stem)
    for x in re.split(r"[|,;\s]+", s):
        x = x.strip()
        if x:
            tokens.add(x)
            tokens.add(Path(x).stem)
    # Common Genome Comparator uploaded-genome labels.
    m = re.search(r"(AZE_[A-Za-z0-9_.-]+)", s)
    if m:
        tokens.add(m.group(1))
        tokens.add(Path(m.group(1)).stem)
    return tokens


def match_reference_labels(
    profile_labels: Iterable[str],
    reference_metadata: pd.DataFrame,
    *,
    reference_id_col: str,
) -> tuple[dict[str, str], pd.DataFrame]:
    """Map profile-matrix labels to reference metadata IDs conservatively."""
    if reference_id_col not in reference_metadata.columns:
        raise ValueError(f"Reference ID column not found: {reference_id_col}")

    labels = [str(x) for x in profile_labels]
    token_map = {lab: _label_tokens(lab) for lab in labels}
    mapping: dict[str, str] = {}
    rows = []

    for ref_id in reference_metadata[reference_id_col].astype(str):
        candidates = [lab for lab, toks in token_map.items() if ref_id in toks]
        if len(candidates) == 0 and ref_id.isdigit():
            pat = re.compile(rf"(?<!\d){re.escape(ref_id)}(?!\d)")
            candidates = [lab for lab in labels if pat.search(lab)]
        status = "matched" if len(candidates) == 1 else ("ambiguous" if len(candidates) > 1 else "unmatched")
        matched = candidates[0] if len(candidates) == 1 else ""
        if matched:
            if matched in mapping:
                raise ValueError(
                    f"Profile label '{matched}' matched more than one reference metadata row."
                )
            mapping[matched] = ref_id
        rows.append(
            {
                reference_id_col: ref_id,
                "profile_label": matched,
                "match_status": status,
                "candidate_labels": ";".join(candidates),
            }
        )

    return mapping, pd.DataFrame(rows)


def _pairwise_distance(
    query: pd.Series,
    reference: pd.Series,
    *,
    missing_tokens: set[str],
) -> tuple[int, int, int, float]:
    """Return raw AD, shared loci, missing-in-either, BIGSdb-style normalised AD."""
    q = np.array([_normalise_allele(v, missing_tokens) for v in query], dtype=object)
    r = np.array([_normalise_allele(v, missing_tokens) for v in reference], dtype=object)
    q_ok = np.array([x is not None for x in q])
    r_ok = np.array([x is not None for x in r])
    common = q_ok & r_ok
    shared = int(common.sum())
    total = len(q)
    if shared == 0:
        return 0, 0, total, float("inf")
    diffs = int((q[common] != r[common]).sum())
    missing_in_either = total - shared
    # BIGSdb lincodes.pl uses 100*diffs/shared_loci and compares this against
    # thresholds converted from allele differences across the full scheme.
    # Converting back to an allele-difference scale gives:
    normalised_ad = float(diffs * total / shared)
    return diffs, shared, missing_in_either, normalised_ad


def _majority(values: Sequence[str]) -> tuple[str, int, float, int]:
    vals = [str(x) for x in values if str(x).strip()]
    if not vals:
        return "", 0, 0.0, 0
    counts = Counter(vals)
    value, n = counts.most_common(1)[0]
    return value, int(n), float(n / len(vals)), len(counts)


def place_profiles(
    profiles: pd.DataFrame,
    reference_metadata: pd.DataFrame,
    *,
    reference_id_col: str = "pubmlst_id",
    lin_col: str = "LINcode_v2",
    cgst_col: str = "cgST_v2",
    query_prefix: str | None = None,
    thresholds: Sequence[int] = CAMPY_V2_LIN_DIFF_THRESHOLDS,
    cgc_thresholds: Sequence[int] = CAMPY_V2_CGC_THRESHOLDS,
    min_support: int = 1,
    min_prop: float = 1.0,
    missing_tokens: set[str] | None = None,
) -> PlacementResult:
    """Place unlabelled profiles relative to LIN-coded reference genomes.

    This does not create official LINcodes. It reports nearest official anchors
    and conservative supported LIN prefixes / cgc2 groups.
    """
    if not 0 < min_prop <= 1:
        raise ValueError("min_prop must be in (0,1].")
    if min_support < 1:
        raise ValueError("min_support must be >=1.")

    thresholds = [int(x) for x in thresholds]
    missing_tokens = {x.lower() for x in (missing_tokens or DEFAULT_MISSING_TOKENS)}

    if lin_col not in reference_metadata.columns:
        raise ValueError(f"Reference LIN column not found: {lin_col}")

    mapping, label_map = match_reference_labels(
        profiles.index,
        reference_metadata,
        reference_id_col=reference_id_col,
    )
    if not mapping:
        raise ValueError("No reference metadata rows could be matched to profile labels.")

    meta = reference_metadata.copy()
    meta[reference_id_col] = meta[reference_id_col].astype(str)
    meta = meta.set_index(reference_id_col, drop=False)

    ref_labels = list(mapping)
    if query_prefix:
        query_labels = [
            x for x in profiles.index
            if str(x).startswith(query_prefix) and x not in mapping
        ]
    else:
        query_labels = [x for x in profiles.index if x not in mapping]
    if not query_labels:
        raise ValueError("No query genomes found after matching reference labels.")

    pairwise_rows = []
    for qid in query_labels:
        q = profiles.loc[qid]
        for rlabel in ref_labels:
            rid = mapping[rlabel]
            raw, shared, missing, norm = _pairwise_distance(
                q,
                profiles.loc[rlabel],
                missing_tokens=missing_tokens,
            )
            pairwise_rows.append(
                {
                    "query_id": str(qid),
                    "reference_matrix_id": str(rlabel),
                    "reference_id": rid,
                    "raw_AD": raw,
                    "shared_loci": shared,
                    "missing_in_either": missing,
                    "normalised_AD": norm,
                }
            )
    pairwise = pd.DataFrame(pairwise_rows)

    threshold_rows = []
    cgc_rows = []
    summary_rows = []

    for qid, d in pairwise.groupby("query_id", sort=False):
        d = d.sort_values(["normalised_AD", "raw_AD", "reference_id"]).reset_index(drop=True)
        min_ad = float(d["normalised_AD"].min())
        nearest = d[np.isclose(d["normalised_AD"], min_ad, rtol=0, atol=1e-9)]
        nearest_ids = nearest["reference_id"].astype(str).tolist()
        nearest_meta = meta.loc[nearest_ids]

        nearest_lins = sorted(
            set(x for x in nearest_meta[lin_col].astype(str) if x and x.lower() != "nan")
        )
        nearest_cgsts = []
        if cgst_col in nearest_meta.columns:
            nearest_cgsts = sorted(
                set(x for x in nearest_meta[cgst_col].astype(str) if x and x.lower() != "nan")
            )

        deepest_level = 0
        deepest_threshold = None
        deepest_prefix = ""
        deepest_support_n = 0
        deepest_prop = 0.0
        ambiguous_finer = False

        for level, threshold in enumerate(thresholds, start=1):
            eligible = d[d["normalised_AD"] <= float(threshold) + 1e-9]
            ids = eligible["reference_id"].astype(str).tolist()
            prefixes = []
            for rid in ids:
                lin = str(meta.loc[rid, lin_col]).strip()
                if not lin or lin.lower() in {"nan", "none", "na"}:
                    continue
                pref = lin_prefix(lin, level)
                if pref:
                    prefixes.append(pref)

            dominant, dom_n, dom_prop, n_unique = _majority(prefixes)
            if not prefixes:
                status = "NO_REFERENCE"
            elif dom_n >= min_support and dom_prop >= min_prop:
                status = "SUPPORTED"
                deepest_level = level
                deepest_threshold = threshold
                deepest_prefix = dominant
                deepest_support_n = dom_n
                deepest_prop = dom_prop
            else:
                status = "AMBIGUOUS"
                if deepest_level:
                    ambiguous_finer = True

            threshold_rows.append(
                {
                    "query_id": qid,
                    "lin_level": level,
                    "difference_threshold": threshold,
                    "eligible_reference_n": len(ids),
                    "lin_labeled_reference_n": len(prefixes),
                    "dominant_prefix": dominant,
                    "dominant_n": dom_n,
                    "dominant_prop": dom_prop,
                    "unique_prefixes": n_unique,
                    "status": status,
                }
            )

        cgc_values = {}
        for threshold in cgc_thresholds:
            col = f"Cjc_cgc2_{threshold}"
            if col not in meta.columns:
                continue
            eligible = d[d["normalised_AD"] <= float(threshold) + 1e-9]
            ids = eligible["reference_id"].astype(str).tolist()
            vals = [
                str(meta.loc[rid, col]).strip()
                for rid in ids
                if str(meta.loc[rid, col]).strip()
                and str(meta.loc[rid, col]).strip().lower() not in {"nan", "none", "na"}
            ]
            dominant, dom_n, dom_prop, n_unique = _majority(vals)
            status = (
                "NO_REFERENCE" if not vals
                else "SUPPORTED" if dom_n >= min_support and dom_prop >= min_prop
                else "AMBIGUOUS"
            )
            cgc_rows.append(
                {
                    "query_id": qid,
                    "difference_threshold": int(threshold),
                    "classification": col,
                    "eligible_reference_n": len(ids),
                    "labeled_reference_n": len(vals),
                    "dominant_group": dominant,
                    "dominant_n": dom_n,
                    "dominant_prop": dom_prop,
                    "unique_groups": n_unique,
                    "status": status,
                }
            )
            cgc_values[f"{col}_placement"] = dominant if status == "SUPPORTED" else ""
            cgc_values[f"{col}_status"] = status
            cgc_values[f"{col}_support_n"] = dom_n
            cgc_values[f"{col}_support_prop"] = dom_prop

        summary_rows.append(
            {
                "query_id": qid,
                "nearest_reference_id": ";".join(nearest_ids),
                "nearest_raw_AD": ";".join(str(x) for x in nearest["raw_AD"].tolist()),
                "nearest_normalised_AD": min_ad,
                "nearest_shared_loci": ";".join(str(x) for x in nearest["shared_loci"].tolist()),
                "nearest_cgST": ";".join(nearest_cgsts),
                "nearest_official_LINcode": ";".join(nearest_lins),
                "deepest_supported_LIN_level": deepest_level if deepest_level else "",
                "deepest_supported_difference_threshold": (
                    deepest_threshold if deepest_threshold is not None else ""
                ),
                "deepest_supported_LIN_prefix": deepest_prefix,
                "LIN_support_n": deepest_support_n,
                "LIN_support_prop": deepest_prop,
                "ambiguous_at_finer_level": "Yes" if ambiguous_finer else "No",
                "placement_status": "SUPPORTED_PREFIX" if deepest_prefix else "UNRESOLVED",
                **cgc_values,
            }
        )

    return PlacementResult(
        summary=pd.DataFrame(summary_rows),
        by_threshold=pd.DataFrame(threshold_rows),
        cgc2=pd.DataFrame(cgc_rows),
        pairwise=pairwise,
        label_map=label_map,
    )


def place_from_files(
    profiles_path: str | Path,
    reference_metadata_path: str | Path,
    *,
    sheet: str = "all",
    sample_id_col: str | None = None,
    locus_prefix: str = "CAMP",
    reference_id_col: str = "pubmlst_id",
    lin_col: str = "LINcode_v2",
    cgst_col: str = "cgST_v2",
    query_prefix: str | None = None,
    thresholds: Sequence[int] = CAMPY_V2_LIN_DIFF_THRESHOLDS,
    min_support: int = 1,
    min_prop: float = 1.0,
) -> PlacementResult:
    profiles = read_profiles(
        profiles_path,
        sheet=sheet,
        sample_id_col=sample_id_col,
        locus_prefix=locus_prefix,
    )
    refs = pd.read_csv(reference_metadata_path, sep="\t", dtype=str, low_memory=False).fillna("")
    return place_profiles(
        profiles,
        refs,
        reference_id_col=reference_id_col,
        lin_col=lin_col,
        cgst_col=cgst_col,
        query_prefix=query_prefix,
        thresholds=thresholds,
        min_support=min_support,
        min_prop=min_prop,
    )
