import pytest
import io
from cite_seq_count import preprocessing


@pytest.fixture
def data():
    from collections import OrderedDict
    from itertools import islice
    
    # Test file paths
    pytest.correct_whitelist_path = 'tests/test_data/whitelists/correct.csv'
    pytest.correct_tags_path = 'tests/test_data/tags/correct.csv'
    pytest.correct_R1_path = 'tests/test_data/fastq/correct_R1.fastq.gz'
    pytest.correct_R2_path = 'tests/test_data/fastq/correct_R2.fastq.gz'
    pytest.corrupt_R1_path = 'tests/test_data/fastq/corrupted_R1.fastq.gz'
    pytest.corrupt_R2_path = 'tests/test_data/fastq/corrupted_R2.fastq.gz'
    
    # Create some variables to compare to
    pytest.correct_whitelist = set(['ACTGTTTTATTGGCCT','TTCATAAGGTAGGGAT'])
    pytest.correct_tags = {
        'AGGACCATCCAA':'CITE_LEN_12_1',
        'ACATGTTACCGT':'CITE_LEN_12_2',
        'AGCTTACTATCC':'CITE_LEN_12_3',
        'TCGATAATGCGAGTACAA':'CITE_LEN_18_1',
        'GAGGCTGAGCTAGCTAGT':'CITE_LEN_18_2',
        'GGCTGATGCTGACTGCTA':'CITE_LEN_18_3',
        'TGTGACGTATTGCTAGCTAG':'CITE_LEN_20_1',
        'ACTGTCTAACGGGTCAGTGC':'CITE_LEN_20_2',
        'TATCACATCGGTGGATCCAT':'CITE_LEN_20_3'}
    pytest.correct_ordered_tags = OrderedDict({
        'TGTGACGTATTGCTAGCTAG':'CITE_LEN_20_1-TGTGACGTATTGCTAGCTAG',
        'ACTGTCTAACGGGTCAGTGC':'CITE_LEN_20_2-ACTGTCTAACGGGTCAGTGC',
        'TATCACATCGGTGGATCCAT':'CITE_LEN_20_3-TATCACATCGGTGGATCCAT',
        'TCGATAATGCGAGTACAA':'CITE_LEN_18_1-TCGATAATGCGAGTACAA',
        'GAGGCTGAGCTAGCTAGT':'CITE_LEN_18_2-GAGGCTGAGCTAGCTAGT',
        'GGCTGATGCTGACTGCTA':'CITE_LEN_18_3-GGCTGATGCTGACTGCTA',
        'AGGACCATCCAA':'CITE_LEN_12_1-AGGACCATCCAA',
        'ACATGTTACCGT':'CITE_LEN_12_2-ACATGTTACCGT',
        'AGCTTACTATCC':'CITE_LEN_12_3-AGCTTACTATCC'})
    pytest.barcode_slice = slice(0, 16)
    pytest.umi_slice = slice(16, 26)
    pytest.barcode_umi_length = 26

@pytest.mark.dependency()
def test_parse_whitelist_csv(data):
    assert preprocessing.parse_whitelist_csv(pytest.correct_whitelist_path, 16, 1) == (pytest.correct_whitelist,1)

@pytest.mark.dependency()
def test_parse_tags_csv(data):
    assert preprocessing.parse_tags_csv(pytest.correct_tags_path) == pytest.correct_tags

@pytest.mark.dependency(depends=['test_parse_tags_csv'])
def test_check_tags(data):
    assert preprocessing.check_tags(pytest.correct_tags, 5) == pytest.correct_ordered_tags

@pytest.mark.dependency(depends=['test_check_tags'])
def test_check_distance_too_big_between_tags(data):
    with pytest.raises(SystemExit):
        preprocessing.check_tags(pytest.correct_tags, 8)

@pytest.mark.dependency(depends=['test_parse_whitelist_csv'])
def test_check_barcodes_lengths(data):
    assert preprocessing.check_barcodes_lengths(26, 1, 16, 17, 26) == (pytest.barcode_slice, pytest.umi_slice, pytest.barcode_umi_length)

@pytest.mark.dependency()
def test_get_n_lines(data):
  assert preprocessing.get_n_lines(pytest.correct_R1_path) == (200 * 4)

@pytest.mark.dependency(depends=['test_get_n_lines'])
def test_get_n_lines_not_multiple_of_4(data):
  with pytest.raises(SystemExit):
    preprocessing.get_n_lines(pytest.corrupt_R1_path)