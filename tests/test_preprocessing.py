"""Test function preprocessing of the module"""

import glob
from collections import namedtuple

import pytest
import polars as pl
from polars.testing import assert_frame_equal
from cite_seq_count import preprocessing


@pytest.fixture
def data():
    """Load up data for testing"""

    pytest.passing_csv = "tests/test_data/tags/pass/*.csv"
    pytest.failing_csv = "tests/test_data/tags/fail/*.csv"

    pytest.passing_reference_list_csv = "tests/test_data/reference_lists/pass/*.csv"
    pytest.failing_reference_list_csv = "tests/test_data/reference_lists/fail/*.csv"

    pytest.passing_filtered_list_csv = "tests/test_data/filtered_lists/pass/*.csv"
    pytest.failing_filtered_list_csv = "tests/test_data/filtered_lists/fail/*.csv"

    pytest.correct_tags_path = "tests/test_data/tags/pass/correct.csv"

    # Create some variables to compare to
    pytest.correct_reference_translation_list = pl.DataFrame(
        {"reference": "ACTGTTTTATTGGCCT", "translation": "TTCATCCTTTAGGGAT"}
    )
    pytest.correct_tag_pl = pl.DataFrame(
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
    pytest.barcode_slice = slice(0, 16)
    pytest.umi_slice = slice(16, 26)
    pytest.barcode_umi_length = 26


def test_csv_parser(data):
    """Test the csv parser

    Args:
        data (_type_): _description_
    """
    passing_files = glob.glob(pytest.passing_csv)
    for file_path in passing_files:
        preprocessing.parse_tags_csv(file_path)
    with pytest.raises(SystemExit):
        failing_files = glob.glob(pytest.failing_csv)
        for file_path in failing_files:
            preprocessing.parse_tags_csv(file_path)


@pytest.mark.dependency()
def test_parse_reference_list_csv(data):
    passing_files = glob.glob(pytest.passing_reference_list_csv)
    for file_path in passing_files:
        assert_frame_equal(
            left=preprocessing.parse_barcode_reference(file_path, 16),
            right=pytest.correct_reference_translation_list,
        )
    with pytest.raises(SystemExit):
        failing_files = glob.glob(pytest.failing_reference_list_csv)
        for file_path in failing_files:
            preprocessing.parse_barcode_reference(file_path, 16)


def test_check_distance_too_big_between_tags(data):
    with pytest.raises(SystemExit):
        preprocessing.check_tags(pytest.correct_tag_pl, 8)
