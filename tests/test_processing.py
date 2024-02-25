import pytest
from cite_seq_count import processing
import polars as pl
import polars_distance as pld
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
        },
        schema={"barcode": pl.String, "count": pl.UInt32},
    )


def test_correct_barcodes_pl(barcodes_df):
    barcode_subset_df = pl.DataFrame(
        {
            "subset": [
                "TACATATTCTTTACTG",
                "GCTAGTCGTAGCTAGA",
                "TAGAGGGAGGTCAAGC",
            ]
        }, schema = {"subset":pl.String}
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
        },
        schema={"barcode": pl.String, "count": pl.UInt32},
    )
    assert_frame_equal(
        corrected_barcodes, expected_corrected_barcodes, check_row_order=False
    )

    # Assert the number of corrected barcodes
    expected_n_corrected_barcodes = 9
    assert n_corrected_barcodes == expected_n_corrected_barcodes

    # Assert the mapped barcodes
    expected_mapped_barcodes = {
        "AACATATTCTTTACTG": "TACATATTCTTTACTG",
        "CACATATTCTTTACTG": "TACATATTCTTTACTG",
        "GACATATTCTTTACTG": "TACATATTCTTTACTG",
        "GCTAGTCGTAGCTAGT": "GCTAGTCGTAGCTAGA",
        "GCTAGTCGTAGCTAGG": "GCTAGTCGTAGCTAGA",
        "GCTAGTCGTAGCTAGC": "GCTAGTCGTAGCTAGA",
        "TAGAGGGACGTCAAGC": "TAGAGGGAGGTCAAGC",
        "TAGAGGGATGTCAAGC": "TAGAGGGAGGTCAAGC",
        "TAGAGGGAAGTCAAGC": "TAGAGGGAGGTCAAGC",
    }
    assert mapped_barcodes == expected_mapped_barcodes


# def test_find_closest_match():
#     input_lf = pl.LazyFrame(
#         {"sequence": ["ATGCTATCAG", "GGCGAGGCT", "GGATTATCGA", "GCTAGCTTAG"]}
#     )
#     ref_seq = pl.LazyFrame(
#         {"ref": ["ATGCTAACAG", "GCTAGCTAT", "AGGAGATC", "GGATAGCGA", "GATTCGGAG"]}
#     )
#     res = processing.find_closest_match(
#         df=input_lf, source_column="sequence", target_df=ref_seq
#     )
#     print(res.collect())
#     assert_frame_equal(res, pl.LazyFrame())


# def test_umi_correction():
#     # Broadest context
#     test_data_broad = pl.DataFrame(
#         {   "barcode": ["cell1", "cell1", "cell1", "cell2", "cell2", "cell2"],
#             "feature": ["gene1", "gene1", "gene2", "gene1", "gene2", "gene2"],
#             "umi": ["CCCC", "CCCA", "TTTT", "AAAA", "GGGG", "CAAA"],
#             "count": [10, 1, 2, 3, 6, 3],
#         }
#     )

#     expected_data_broad = pl.DataFrame(
#         {   "barcode": ["cell1", "cell1", "cell2", "cell2"],
#             "feature": ["gene1", "gene2", "gene1", "gene2"],
#             "umi": ["CCCC", "TTTT", "AAAA", "GGGG"],
#             "count": [11, 2, 6, 6],
#         }
#     )
#     # One cell, one feature
#     test_data = pl.DataFrame(
#         {
#             "umi": ["CCCC", "CCCA", "TTTT"],
#             "count": [10, 1, 2],
#         }
#     )
#     expected = pl.DataFrame(
#         {
#             "umi": ["CCCC", "TTTT"],
#             "count": [11, 2],
#         }
#     )
#     # What I've got so far that deals with one case of the smaller context but doesn't keep the "uncorrected" umis
#     res = (
#         test_data.join(test_data.select("umi"), on="umi", how="cross")
#         .filter(pl.col("umi") != pl.col("umi_right"))
#         .with_columns(pld.col("umi").dist_str.hamming("umi_right").alias("hamming"))
#         .filter(pl.col("hamming") <= 1)
#         .sort("count", descending=True)
#         .select(pl.first("umi"), pl.sum("count"))
#     )
#     print(res)
