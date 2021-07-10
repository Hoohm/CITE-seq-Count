import pytest
import io
from cite_seq_count import preprocessing
import glob
from collections import namedtuple, OrderedDict
from itertools import islice


@pytest.fixture
def data():

    pytest.passing_csv = "tests/test_data/tags/pass/*.csv"
    pytest.failing_csv = "tests/test_data/tags/fail/*.csv"

    pytest.passing_reference_list_csv = "tests/test_data/reference_lists/pass/*.csv"
    pytest.failing_reference_list_csv = "tests/test_data/reference_lists/fail/*.csv"

    pytest.passing_filtered_list_csv = "tests/test_data/filtered_lists/pass/*.csv"
    pytest.failing_filtered_list_csv = "tests/test_data/filtered_lists/fail/*.csv"

    pytest.correct_tags_path = "tests/test_data/tags/pass/correct.csv"

    # Create some variables to compare to
    pytest.correct_reference_list = set(["ACTGTTTTATTGGCCT", "TTCATAAGGTAGGGAT"])
    pytest.correct_tags = {
        "AGGACCATCCAA": "CITE_LEN_12_1",
        "ACATGTTACCGT": "CITE_LEN_12_2",
        "AGCTTACTATCC": "CITE_LEN_12_3",
        "TCGATAATGCGAGTACAA": "CITE_LEN_18_1",
        "GAGGCTGAGCTAGCTAGT": "CITE_LEN_18_2",
        "GGCTGATGCTGACTGCTA": "CITE_LEN_18_3",
        "TGTGACGTATTGCTAGCTAG": "CITE_LEN_20_1",
        "ACTGTCTAACGGGTCAGTGC": "CITE_LEN_20_2",
        "TATCACATCGGTGGATCCAT": "CITE_LEN_20_3",
    }
    tag = namedtuple("tag", ["name", "sequence", "id"])
    pytest.correct_tags_tuple = [
        tag(name="CITE_LEN_20_1", sequence="TGTGACGTATTGCTAGCTAG", id=0),
        tag(name="CITE_LEN_20_2", sequence="ACTGTCTAACGGGTCAGTGC", id=1),
        tag(name="CITE_LEN_20_3", sequence="TATCACATCGGTGGATCCAT", id=2),
        tag(name="CITE_LEN_18_1", sequence="TCGATAATGCGAGTACAA", id=3),
        tag(name="CITE_LEN_18_2", sequence="GAGGCTGAGCTAGCTAGT", id=4),
        tag(name="CITE_LEN_18_3", sequence="GGCTGATGCTGACTGCTA", id=5),
        tag(name="CITE_LEN_12_1", sequence="AGGACCATCCAA", id=6),
        tag(name="CITE_LEN_12_2", sequence="ACATGTTACCGT", id=7),
        tag(name="CITE_LEN_12_3", sequence="AGCTTACTATCC", id=8),
    ]
    pytest.barcode_slice = slice(0, 16)
    pytest.umi_slice = slice(16, 26)
    pytest.barcode_umi_length = 26


def test_csv_parser(data):
    passing_files = glob.glob(pytest.passing_csv)
    for file_path in passing_files:
        preprocessing.parse_tags_csv(file_path)
    with pytest.raises(SystemExit):
        failing_files = glob.glob(pytest.failing_csv)
        for file_path in failing_files:
            preprocessing.parse_tags_csv(file_path)


def test_filtered_list_parser(data):
    passing_files = glob.glob(pytest.passing_filtered_list_csv)
    for file_path in passing_files:
        preprocessing.parse_filtered_list_csv(file_path, barcode_length=16)
    with pytest.raises(SystemExit):
        failing_files = glob.glob(pytest.failing_filtered_list_csv)
        for file_path in failing_files:
            preprocessing.parse_filtered_list_csv(file_path, barcode_length=16)


@pytest.mark.dependency()
def test_parse_reference_list_csv(data):
    passing_files = glob.glob(pytest.passing_reference_list_csv)
    for file_path in passing_files:
        assert preprocessing.parse_cell_list_csv(file_path, 16).keys() in (
            pytest.correct_reference_list,
            1,
        )
    with pytest.raises(SystemExit):
        failing_files = glob.glob(pytest.failing_reference_list_csv)
        for file_path in failing_files:
            preprocessing.parse_cell_list_csv(file_path, 16)


@pytest.mark.dependency()
def test_parse_tags_csv(data):
    tags = preprocessing.check_tags(pytest.correct_tags, 5)[0]
    for i, tag in enumerate(tags):
        assert tag == pytest.correct_tags_tuple[i]


@pytest.mark.dependency(depends=["test_parse_tags_csv"])
def test_check_distance_too_big_between_tags(data):
    with pytest.raises(SystemExit):
        preprocessing.check_tags(pytest.correct_tags, 8)

