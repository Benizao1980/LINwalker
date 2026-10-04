import pandas as pd

from linwalker.place import place_profiles
from linwalker.utils import infer_lin_levels


def test_infer_current_campylobacter_v2_depth():
    s = pd.Series([
        "0_2_0_1_0_0_1_0_0_2_0_0_1_0_1_0_0_0",
        None,
    ])
    assert infer_lin_levels(s) == 18


def test_reference_anchored_prefix_placement():
    profiles = pd.DataFrame(
        {
            "CAMP0001": ["1", "2", "1"],
            "CAMP0002": ["1", "2", "1"],
            "CAMP0003": ["1", "2", "1"],
            "CAMP0004": ["1", "2", "2"],
        },
        index=["100", "200", "AZE_Q1"],
    )

    refs = pd.DataFrame(
        {
            "pubmlst_id": ["100", "200"],
            "LINcode_v2": ["0_1_0_0", "0_2_0_0"],
            "cgST_v2": ["10", "20"],
            "Cjc_cgc2_5": ["A", "B"],
        }
    )

    result = place_profiles(
        profiles,
        refs,
        thresholds=[3, 2, 1, 0],
        query_prefix="AZE_",
    )

    row = result.summary.iloc[0]
    assert row["query_id"] == "AZE_Q1"
    assert row["nearest_reference_id"] == "100"
    assert float(row["nearest_normalised_AD"]) == 1.0
    assert row["deepest_supported_difference_threshold"] == 1
    assert row["deepest_supported_LIN_prefix"] == "0_1_0"
    assert row["nearest_cgST"] == "10"
    assert row["placement_status"] == "SUPPORTED_PREFIX"


def test_missing_loci_are_normalised_like_bigsdb():
    profiles = pd.DataFrame(
        {
            "CAMP0001": ["1", "1"],
            "CAMP0002": ["1", "2"],
            "CAMP0003": ["1", "X"],
            "CAMP0004": ["1", "1"],
        },
        index=["100", "AZE_Q1"],
    )
    refs = pd.DataFrame(
        {
            "pubmlst_id": ["100"],
            "LINcode_v2": ["0_0_0_0"],
            "cgST_v2": ["10"],
        }
    )

    result = place_profiles(
        profiles,
        refs,
        thresholds=[3, 2, 1, 0],
        query_prefix="AZE_",
    )

    pair = result.pairwise.iloc[0]
    assert pair["raw_AD"] == 1
    assert pair["shared_loci"] == 3
    assert abs(pair["normalised_AD"] - (4 / 3)) < 1e-9


def test_exact_reference_profile_reports_official_cgst_and_context():
    profiles = pd.DataFrame(
        {
            "CAMP0001": ["1", "1"],
            "CAMP0002": ["2", "2"],
            "CAMP0003": ["3", "3"],
        },
        index=["100", "AZE_EXACT"],
    )
    refs = pd.DataFrame(
        {
            "pubmlst_id": ["100"],
            "LINcode_v2": ["0_2_0"],
            "cgST_v2": ["72756"],
            "ST": ["10042"],
            "clonal_complex": ["ST-828 complex"],
        }
    )

    result = place_profiles(
        profiles,
        refs,
        thresholds=[2, 1, 0],
        query_prefix="AZE_",
    )

    row = result.summary.iloc[0]
    assert row["placement_status"] == "EXACT_REFERENCE_PROFILE"
    assert row["exact_official_cgST"] == "72756"
    assert row["nearest_ST"] == "10042"
    assert row["nearest_clonal_complex"] == "ST-828 complex"
    assert float(row["nearest_normalised_AD"]) == 0.0
