import pytest
from cite_seq_count import processing
import polars as pl
from polars.testing import assert_frame_equal


@pytest.fixture
def barcodes_df():
    return pl.DataFrame(
        {
            "barcode": [
                "TACATATTCTTTACTG",
                "AACATATTCTTTACTG",
                "CACATATTCTTTACTG",
                "GACATATTCTTTACTG",
                "GCTAGTCGTAGCTAGA",
                "GCTAGTCGTAGCTAGT",
                "GCTAGTCGTAGCTAGG",
                "GCTAGTCGTAGCTAGC",
                "TAGAGGGAGGTCAAGC",
                "TAGAGGGACGTCAAGC",
                "TAGAGGGATGTCAAGC",
                "TAGAGGGAAGTCAAGC",
            ],
            "count": [5, 1, 1, 1, 5, 1, 1, 1, 5, 1, 1, 1],
        }
    )


def test_correct_barcodes_pl(barcodes_df):
    barcode_subset_df = pl.DataFrame(
        {
            "subset": [
                "TACATATTCTTTACTG",
                "GCTAGTCGTAGCTAGA",
                "TAGAGGGAGGTCAAGC",
            ],
        }
    )
    hamming_distance = 1

    (
        corrected_barcodes,
        n_corrected_barcodes,
        mapped_barcodes,
    ) = processing.correct_barcodes_pl(barcodes_df, barcode_subset_df, hamming_distance)

    # Assert the corrected barcodes
    expected_corrected_barcodes = pl.LazyFrame(
        {
            "barcode": [
                "TACATATTCTTTACTG",
                "GCTAGTCGTAGCTAGA",
                "TAGAGGGAGGTCAAGC",
            ],
            "count": [8, 8, 8],
        }
    )
    assert_frame_equal(
        corrected_barcodes, expected_corrected_barcodes, check_row_order=False
    )

    # Assert the number of corrected barcodes
    expected_n_corrected_barcodes = 9
    assert n_corrected_barcodes == expected_n_corrected_barcodes

    # Assert the mapped barcodes
    expected_mapped_barcodes = {
        "AACATATTCTTTACTG":"TACATATTCTTTACTG",
        "CACATATTCTTTACTG":"TACATATTCTTTACTG",
        "GACATATTCTTTACTG":"TACATATTCTTTACTG",
        "GCTAGTCGTAGCTAGT":"GCTAGTCGTAGCTAGA",
        "GCTAGTCGTAGCTAGG":"GCTAGTCGTAGCTAGA",
        "GCTAGTCGTAGCTAGC":"GCTAGTCGTAGCTAGA",
        "TAGAGGGACGTCAAGC":"TAGAGGGAGGTCAAGC",
        "TAGAGGGATGTCAAGC":"TAGAGGGAGGTCAAGC",
        "TAGAGGGAAGTCAAGC":"TAGAGGGAGGTCAAGC",
    }
    assert mapped_barcodes == expected_mapped_barcodes
