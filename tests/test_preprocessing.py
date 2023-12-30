"""Test function preprocessing of the module"""

import glob

import pytest
import polars as pl
from polars.testing import assert_frame_equal
from cite_seq_count.preprocessing import (
    parse_barcode_file,
    parse_tags_csv,
    check_tags,
    find_knee_estimated_barcodes,
    get_barcode_subset,
)
from cite_seq_count.constants import (
    BARCODE_COLUMN,
    COUNT_COLUMN,
    REFERENCE_COLUMN,
    SUBSET_COLUMN,
    TRANSLATION_COLUMN,
)
from cite_seq_count.chemistry import Chemistry


@pytest.fixture
def passing_references():
    return "tests/test_data/reference_lists/pass/*.csv"


@pytest.fixture
def failing_references():
    return "tests/test_data/reference_lists/fail/*.csv"


@pytest.fixture
def passing_subsets():
    return "tests/test_data/reference_subsets/pass/*.csv"


@pytest.fixture
def failing_subsets():
    return "tests/test_data/reference_subsets/fail/*.csv"


@pytest.fixture
def passing_tags():
    return "tests/test_data/tags/pass/*.csv"


@pytest.fixture
def failing_tags():
    return "tests/test_data/tags/fail/*.csv"


@pytest.fixture
def correct_reference_df():
    return pl.DataFrame(
        {
            REFERENCE_COLUMN: ["ATGCCC", "ATGCTT", "CCGCCC"],
            TRANSLATION_COLUMN: ["GGATCG", "GGATCA", "GATCAA"],
        }
    )


@pytest.fixture
def correct_subset_df():
    return pl.DataFrame({REFERENCE_COLUMN: ["ATGCCC", "ATGCTT"]})


@pytest.fixture
def barcodes_df():
    return pl.DataFrame(
        {
            BARCODE_COLUMN: ["ATGCCC", "ATGCTT", "CCGCCC", "ATATCC", "ATATGG"],
            COUNT_COLUMN: [200, 200, 200, 20, 10],
        }
    )


@pytest.fixture
def correct_tags_df():
    return pl.DataFrame(
        {
            "feature_name": [
                "CITE_LEN_20_1",
                "CITE_LEN_20_2",
                "CITE_LEN_20_3",
                "CITE_LEN_18_1",
                "CITE_LEN_18_2",
                "CITE_LEN_18_3",
                "CITE_LEN_12_1",
                "CITE_LEN_12_2",
                "CITE_LEN_12_3",
            ],
            "sequence": [
                "TGTGACGTATTGCTAGCTAG",
                "ACTGTCTAACGGGTCAGTGC",
                "TATCACATCGGTGGATCCAT",
                "TCGATAATGCGAGTACAA",
                "GAGGCTGAGCTAGCTAGT",
                "GGCTGATGCTGACTGCTA",
                "AGGACCATCCAA",
                "ACATGTTACCGT",
                "AGCTTACTATCC",
            ],
        }
    )


@pytest.fixture
def chemistry_def():
    return Chemistry(
        name="test",
        cell_barcode_start=1,
        cell_barcode_end=6,
        umi_barcode_start=7,
        umi_barcode_end=12,
        r2_trim_start=0,
        barcode_reference_path="tests/test_data/reference_lists/pass/translation.csv",
    )


def test_passing_parse_reference_list_csv(passing_references, correct_reference_df):
    passing_files = glob.glob(passing_references)
    for file_path in passing_files:
        assert_frame_equal(
            left=parse_barcode_file(file_path, 16, [REFERENCE_COLUMN]),
            right=correct_reference_df,
        )


def test_failing_parse_reference_list_csv(failing_references):
    with pytest.raises(SystemExit):
        failing_files = glob.glob(failing_references)
        for file_path in failing_files:
            parse_barcode_file(file_path, 16, [REFERENCE_COLUMN])


def test_parse_subset_list_csv(passing_subsets, failing_subsets, correct_subset_df):
    passing_files = glob.glob(passing_subsets)
    for file_path in passing_files:
        assert_frame_equal(
            left=parse_barcode_file(file_path, 16, [REFERENCE_COLUMN]),
            right=correct_subset_df,
        )
    with pytest.raises(SystemExit):
        failing_files = glob.glob(failing_subsets)
        for file_path in failing_files:
            parse_barcode_file(file_path, 16, [REFERENCE_COLUMN])


def test_check_distance_too_big_between_tags(correct_tags_df):
    with pytest.raises(SystemExit):
        check_tags(correct_tags_df, 8)


# Test if there is no reference and no whitelist
def test_find_knee_estimated_barcodes(barcodes_df):
    expected_subset = pl.DataFrame({SUBSET_COLUMN: ["ATGCCC", "ATGCTT", "CCGCCC"]})
    subset = find_knee_estimated_barcodes(barcodes_df=barcodes_df)
    assert_frame_equal(subset, expected_subset)


@pytest.fixture
def barcode_subset_df():
    return pl.DataFrame({SUBSET_COLUMN: ["ATGCCC", "ATGCTT"]})


def test_get_barcode_subset_with_reference(correct_reference_df, barcodes_df):
    expected_subset = pl.DataFrame(
        {
            SUBSET_COLUMN: ["ATGCCC", "ATGCTT"],
        }
    )
    expected_enable_correction = True

    subset, enable_correction = get_barcode_subset(
        barcode_reference=correct_reference_df,
        n_barcodes=2,
        chemistry=None,
        barcode_subset=None,
        barcodes_df=barcodes_df,
    )

    assert_frame_equal(subset, expected_subset)
    assert enable_correction == expected_enable_correction


def test_get_barcode_subset_without_reference(barcodes_df):
    expected_subset = pl.DataFrame({SUBSET_COLUMN: ["ATGCCC", "ATGCTT", "CCGCCC"]})
    expected_enable_correction = True

    subset, enable_correction = get_barcode_subset(
        barcode_reference=None,
        n_barcodes=3,
        chemistry=None,
        barcode_subset=None,
        barcodes_df=barcodes_df,
    )

    assert_frame_equal(subset, expected_subset)
    assert enable_correction == expected_enable_correction


def test_get_barcode_subset_with_large_n_barcodes(correct_reference_df, barcodes_df):
    expected_subset = pl.DataFrame(
        {
            SUBSET_COLUMN: ["ATGCCC", "ATGCTT", "CCGCCC"],
        }
    )
    expected_enable_correction = False

    subset, enable_correction = get_barcode_subset(
        barcode_reference=correct_reference_df,
        n_barcodes=4,
        chemistry=None,
        barcode_subset=None,
        barcodes_df=barcodes_df,
    )

    assert_frame_equal(subset, expected_subset)
    assert enable_correction == expected_enable_correction


def test_get_barcode_subset_with_existing_subset(
    correct_reference_df, barcode_subset_df, barcodes_df
):
    expected_subset = pl.DataFrame({SUBSET_COLUMN: ["ATGCCC", "ATGCTT"]})
    expected_enable_correction = True

    subset, enable_correction = get_barcode_subset(
        barcode_reference=correct_reference_df,
        n_barcodes=2,
        chemistry=None,
        barcode_subset=barcode_subset_df,
        barcodes_df=barcodes_df,
    )

    assert_frame_equal(subset, expected_subset)
    assert enable_correction == expected_enable_correction


def test_get_barcode_subset_with_no_subset(correct_reference_df, barcodes_df):
    expected_subset = pl.DataFrame(
        {
            SUBSET_COLUMN: ["ATGCCC", "ATGCTT", "CCGCCC"],
        }
    )
    expected_enable_correction = True

    subset, enable_correction = get_barcode_subset(
        barcode_reference=correct_reference_df,
        n_barcodes=3,
        chemistry=None,
        barcode_subset=None,
        barcodes_df=barcodes_df,
    )

    assert_frame_equal(subset, expected_subset)
    assert enable_correction == expected_enable_correction


def test_parse_tags_csv(correct_tags_df):
    file_name = "tests/test_data/tags/pass/correct.csv"
    result_df = parse_tags_csv(file_name)

    assert_frame_equal(
        result_df, correct_tags_df, check_column_order=False, check_row_order=False
    )
