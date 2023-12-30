"""Mapping module. Holds all code related to mapping reads
"""
import polars as pl
from rapidfuzz import fuzz, process, distance

from cite_seq_count.constants import (
    COUNT_COLUMN,
    SEQUENCE_COLUMN,
    R2_COLUMN,
    FEATURE_NAME_COLUMN,
    UNMAPPED_NAME,
)


def find_best_match_fast(tag_seq, tags_list, maximum_distance):
    choices = tags_list[SEQUENCE_COLUMN].to_list()
    features = tags_list[FEATURE_NAME_COLUMN].to_list()
    res = process.extractOne(choices=choices, query=tag_seq, scorer=fuzz.QRatio)
    min_score = (len(tag_seq) - maximum_distance) / len(tag_seq) * 100
    if res[1] >= min_score:
        return features[res[2]]
    return UNMAPPED_NAME


def map_reads_hybrid(
    r2_df: pl.DataFrame, parsed_tags: pl.DataFrame, maximum_distance: int
) -> tuple[pl.DataFrame, pl.DataFrame]:
    """Map sequence data to a tags reference.
    Using a hybdrid approach where we first join all the data for the exact matches
    then using a hamming distance calculation to find the closest match

    Args:
        r2_df (pl.DataFrame): All r2 sequences to map
        parsed_tags (pl.DataFrame): tags to map to
        maximum_distance (int): max distance allowed for mismatches

    Returns:
        pl.DataFrame: Mapped data
    """
    print("Mapping reads")
    mapped_r2_df = r2_df.join(
        parsed_tags, left_on=R2_COLUMN, right_on=SEQUENCE_COLUMN, how="left"
    ).with_columns(
        pl.when(pl.col(FEATURE_NAME_COLUMN).is_null())
        .then(
            pl.col(R2_COLUMN)
            .map_elements(
                lambda x: find_best_match_fast(
                    x, tags_list=parsed_tags, maximum_distance=maximum_distance
                )
            )
            .alias(FEATURE_NAME_COLUMN)
        )
        .otherwise(pl.col(FEATURE_NAME_COLUMN))
    )
    unmapped_r2_df = mapped_r2_df.filter(pl.col(FEATURE_NAME_COLUMN) == UNMAPPED_NAME)
    print("Mapping done")
    return mapped_r2_df, unmapped_r2_df


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
