"""Mapping module. Holds all code related to mapping reads
"""
from pathlib import Path
import polars as pl
from rapidfuzz import fuzz, process, distance

from cite_seq_count.preprocessing import (
    SEQUENCE_COLUMN,
    R2_COLUMN,
    FEATURE_NAME_COLUMN,
    BARCODE_COLUMN,
    UMI_COLUMN,
    UNMAPPED_NAME,
)


def find_best_match_rapid(tag_seq, tags_list, maximum_distance):
    choices = tags_list[SEQUENCE_COLUMN].to_list()
    features = tags_list[FEATURE_NAME_COLUMN].to_list()
    res = process.extractOne(choices=choices, query=tag_seq, scorer=fuzz.QRatio)
    min_score = (len(tag_seq) - maximum_distance) / len(tag_seq) * 100
    if res[1] >= min_score:
        return features[res[2]]
    return UNMAPPED_NAME


def match_generic_string_dfs(
    ref_df: pl.DataFrame,
    target_df: pl.DataFrame,
    left_on: str,
    right_on: str,
    hamming_distance: int,
):
    corrected_column = "corrected_" + left_on
    joined = (
        target_df.sort(left_on)
        .join_asof(ref_df.sort(right_on), left_on=left_on, right_on=right_on)
        .with_columns(
            pl.when(pl.col(left_on) == pl.col(right_on))
            .then(False)
            .otherwise(True)
            .alias(corrected_column)
        )
        .with_columns(
            pl.when(pl.col(corrected_column))
            .then(distance.Hamming.distance(s1=pl.col(left_on), s2=pl.col(right_on)))
            .otherwise(0)
            .alias("hamming_distance")
        )
    )
    return joined


def map_reads_hybrid(
    r2_df: pl.DataFrame, parsed_tags: pl.DataFrame, maximum_distance: int
) -> pl.DataFrame:
    print("Mapping reads")
    mapped_r2_df = r2_df.join(
        parsed_tags, left_on=R2_COLUMN, right_on=SEQUENCE_COLUMN, how="left"
    ).with_columns(
        pl.when(pl.col(FEATURE_NAME_COLUMN).is_null())
        .then(
            pl.col(R2_COLUMN)
            .map_elements(
                lambda x: find_best_match_rapid(
                    x, tags_list=parsed_tags, maximum_distance=maximum_distance
                )
            )
            .alias(FEATURE_NAME_COLUMN)
        )
        .otherwise(pl.col(FEATURE_NAME_COLUMN))
    )
    print("Mapping done")
    return mapped_r2_df


def check_unmapped(mapped_reads: pl.DataFrame):
    """_summary_

    Args:
        mapped_reads (pl.DataFrame): _description_

    """
    n_reads = mapped_reads.shape[0]
    n_unmapped = mapped_reads.filter(
        pl.col(FEATURE_NAME_COLUMN) == UNMAPPED_NAME
    ).shape[0]
    if n_reads / n_unmapped > 0.99:
        SystemExit("Number of unmapped reads is more than 99%. Exiting")
