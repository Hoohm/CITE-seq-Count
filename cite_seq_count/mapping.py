"""Mapping module. Holds all code related to mapping reads
"""
import polars as pl
from rapidfuzz import fuzz, process
import polars_distance as pld

from cite_seq_count.constants import (
    COUNT_COLUMN,
    SEQUENCE_COLUMN,
    R2_COLUMN,
    FEATURE_NAME_COLUMN,
    UNMAPPED_NAME,
)


def find_best_match_fast(tag_seq, tags_df, maximum_distance):
    min_scores = (
        tags_df.with_columns(pl.col("sequence").str.len_chars().alias("seq_length"))
        .with_columns(
            (
                (pl.col("seq_length") - maximum_distance) / pl.col("seq_length") * 100
            ).alias("min_score")
        )["min_score"]
        .to_list()
    )
    choices = tags_df[SEQUENCE_COLUMN].to_list()
    features = tags_df[FEATURE_NAME_COLUMN].to_list()
    _, score, index = process.extractOne(
        choices=choices, query=tag_seq, scorer=fuzz.partial_ratio
    )
    if score >= min_scores[index]:
        return features[index]
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
                    x, tags_df=parsed_tags, maximum_distance=maximum_distance
                )
            )
            .alias(FEATURE_NAME_COLUMN)
        )
        .otherwise(pl.col(FEATURE_NAME_COLUMN))
    )
    unmapped_r2_df = mapped_r2_df.filter(pl.col(FEATURE_NAME_COLUMN) == UNMAPPED_NAME)
    print("Mapping done")
    return mapped_r2_df, unmapped_r2_df


def map_reads_polars(
    r2_df: pl.DataFrame, parsed_tags: pl.DataFrame, maximum_distance: int
) -> tuple[pl.DataFrame, pl.DataFrame]:
    joined = r2_df.join(
        parsed_tags, left_on=R2_COLUMN, right_on=SEQUENCE_COLUMN, how="left"
    )

    simple_join = joined.filter(~pl.col(FEATURE_NAME_COLUMN).is_null())
    levenshtein_mapped = (
        joined.filter(pl.col(FEATURE_NAME_COLUMN).is_null())
        .join(parsed_tags, left_on=R2_COLUMN, right_on=SEQUENCE_COLUMN, how="cross")
        .drop(FEATURE_NAME_COLUMN)
        .with_columns(
            pld.col(SEQUENCE_COLUMN)
            .dist_str.levenshtein(pl.col(R2_COLUMN))
            .alias("levenshtein_dist")
        )
        .filter(pl.col("levenshtein_dist") <= maximum_distance)
        .drop([SEQUENCE_COLUMN, "levenshtein_dist"])
        .rename({"feature_name_right": FEATURE_NAME_COLUMN})
    )
    multi_mapped = (
        levenshtein_mapped.group_by(R2_COLUMN)
        .agg(pl.count())
        .filter(pl.col(COUNT_COLUMN) > 1)
    )
    mapped = pl.concat(
        [
            simple_join,
            levenshtein_mapped.filter(
                ~pl.col(R2_COLUMN).is_in(multi_mapped[R2_COLUMN])
            ),
        ]
    )
    unmapped = joined.filter(~pl.col(R2_COLUMN).is_in(mapped[R2_COLUMN])).with_columns(
        pl.col(FEATURE_NAME_COLUMN).fill_null(UNMAPPED_NAME)
    )
    return mapped, unmapped


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
