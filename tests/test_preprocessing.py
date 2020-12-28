import pytest
import io
from cite_seq_count import preprocessing
import glob
from collections import namedtuple


@pytest.fixture
def data():
    from collections import OrderedDict
    from itertools import islice

    pytest.passing_csv = "tests/test_data/tags/pass/*.csv"
    pytest.failing_csv = "tests/test_data/tags/fail/*.csv"

    # Test file paths
    pytest.correct_reference_list_path = "tests/test_data/reference_lists/correct.csv"
    pytest.correct_tags_path = "tests/test_data/tags/pass/correct.csv"
    pytest.correct_R1_path = "tests/test_data/fastq/correct_R1.fastq.gz"
    pytest.correct_R2_path = "tests/test_data/fastq/correct_R2.fastq.gz"
    pytest.corrupt_R1_path = "tests/test_data/fastq/corrupted_R1.fastq.gz"
    pytest.corrupt_R2_path = "tests/test_data/fastq/corrupted_R2.fastq.gz"

    pytest.correct_R1_multipath = "path/to/R1_1.fastq.gz,path/to/R1_2.fastq.gz"
    pytest.correct_R2_multipath = "path/to/R2_1.fastq.gz,path/to/R2_2.fastq.gz"
    pytest.incorrect_R2_multipath = (
        "path/to/R2_1.fastq.gz,path/to/R2_2.fastq.gz,path/to/R2_3.fastq.gz"
    )

    pytest.correct_multipath_result = (
        ["path/to/R1_1.fastq.gz", "path/to/R1_2.fastq.gz"],
        ["path/to/R2_1.fastq.gz", "path/to/R2_2.fastq.gz"],
    )

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
        tag(name="unmapped", sequence="UNKNOWN", id=9),
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
            print(file_path)
            preprocessing.parse_tags_csv(file_path)


@pytest.mark.dependency()
def test_parse_reference_list_csv(data):
    assert preprocessing.parse_reference_list_csv(
        pytest.correct_reference_list_path, 16
    ).keys() in (pytest.correct_reference_list, 1,)


@pytest.mark.dependency()
def test_parse_tags_csv(data):
    tags = preprocessing.check_tags(pytest.correct_tags, 5)[0]
    for i, tag in enumerate(tags):
        assert tag == pytest.correct_tags_tuple[i]


@pytest.mark.dependency(depends=["test_parse_tags_csv"])
def test_check_distance_too_big_between_tags(data):
    with pytest.raises(SystemExit):
        preprocessing.check_tags(pytest.correct_tags, 8)


@pytest.mark.dependency()
def test_get_n_lines(data):
    assert preprocessing.get_n_lines(pytest.correct_R1_path) == (200 * 4)


@pytest.mark.dependency(depends=["test_get_n_lines"])
def test_get_n_lines_not_multiple_of_4(data):
    with pytest.raises(SystemExit):
        preprocessing.get_n_lines(pytest.corrupt_R1_path)


@pytest.mark.dependency()
def test_corrrect_multipath(data):
    assert (
        preprocessing.get_read_paths(
            pytest.correct_R1_multipath, pytest.correct_R2_multipath
        )
        == pytest.correct_multipath_result
    )


@pytest.mark.dependency(depends=["test_get_n_lines"])
def test_incorrrect_multipath(data):
    with pytest.raises(SystemExit):
        preprocessing.get_read_paths(
            pytest.correct_R1_multipath, pytest.incorrect_R2_multipath
        )
