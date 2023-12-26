import pytest
from collections import namedtuple
from cite_seq_count import processing
import polars as pl
from polars.testing import assert_frame_equal


@pytest.fixture
def data():
    tag = namedtuple("tag", ["name", "sequence", "id"])

    pytest.barcodes_df = pl.DataFrame(
        {
            "barcode": [
                "TACATATTCTTTACTG",
                "AACATATTCTTTACTG",
                "CACATATTCTTTACTG",
                "GACATATTCTTTACTG",
                "TACATATTCTTTACTA",
                "TACATATTCTTTACTC",
                "TACATATTCTTTACTT",
                "TAGAGGGAGGTCAAGC",
                "TAGAGGGACGTCAAGC",
                "TAGAGGGATGTCAAGC",
                "TAGAGGGAAGTCAAGC",
            ],
            "count": [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
        }
    )

    pytest.barcode_subset_df = pl.DataFrame(
        {"whitelist": ["TACATATTCTTTACTG", "TAGAGGGAAGTCAAGC"]}
    )

    pytest.corrected_barcodes_df = pl.DataFrame(
        {
            "barcode": [
                "TAGAGGGAAGTCAAGC",
                "TACATATTCTTTACTG",
            ],
            "count": [4, 7],
        }
    )

    pytest.expected_cells = 2
    pytest.collapsing_threshold = 1
    pytest.max_umis = 20000


@pytest.mark.dependency()
def test_correct_barcodes(data):
    corrected_barcodes, _, _ = processing.correct_barcodes_pl(
        barcodes_df=pytest.barcodes_df,
        barcode_subset_df=pytest.barcode_subset_df,
        hamming_distance=1,
    )
    assert_frame_equal(
        pytest.corrected_barcodes_df, corrected_barcodes, check_row_order=False
    )


@pytest.mark.dependency()
def test_correct_umis(data):
    temp = processing.correct_umis_in_cells((pytest.results, 2, pytest.max_umis, 2))
    results = temp[0]
    n_corrected = temp[1]
    for cell_barcode in results.keys():
        for TAG in results[cell_barcode]:
            assert len(results[cell_barcode][TAG]) == len(
                pytest.corrected_results[cell_barcode][TAG]
            )
            assert sum(results[cell_barcode][TAG].values()) == sum(
                pytest.corrected_results[cell_barcode][TAG].values()
            )
    assert n_corrected == 3


@pytest.mark.dependency(depends=["test_correct_umis"])
def test_generate_sparse_umi_matrices(data):
    umi_results_matrix = processing.generate_sparse_matrices(
        pytest.corrected_results,
        pytest.tags,
        ["ACTGTTTTATTGGCCT", "TTCATAAGGTAGGGAT"],
        umi_counts=True,
    )
    assert umi_results_matrix.shape == (2, 2)
    total_umis = 0
    for i in range(umi_results_matrix.shape[0]):
        for j in range(umi_results_matrix.shape[1]):
            total_umis += umi_results_matrix[i, j]
    assert total_umis == 3


@pytest.mark.dependency(depends=["test_correct_umis"])
def test_generate_sparse_read_matrices(data):
    read_results_matrix = processing.generate_sparse_matrices(
        pytest.corrected_results,
        pytest.tags,
        ["ACTGTTTTATTGGCCT", "TTCATAAGGTAGGGAT"],
        umi_counts=False,
    )
    assert read_results_matrix.shape == (3, 2)
    total_umis = 0
    for i in range(read_results_matrix.shape[0]):
        for j in range(read_results_matrix.shape[1]):
            total_umis += read_results_matrix[i, j]
    assert total_umis == 12
