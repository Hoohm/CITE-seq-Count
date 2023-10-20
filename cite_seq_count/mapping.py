"""Mapping module. Holds all code related to mapping reads
"""
from pathlib import Path
import polars as pl
from rapidfuzz import fuzz, process

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


def map_reads_hybrid(
    mapping_input_file: Path, parsed_tags: pl.DataFrame, maximum_distance: int
) -> pl.DataFrame:
    input_reads = pl.read_csv(
        mapping_input_file,
        has_header=False,
        new_columns=[BARCODE_COLUMN, UMI_COLUMN, R2_COLUMN],
    )
    mapped_reads = input_reads.join(
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
    return mapped_reads


def map_reads_pl_beta(
    parsed_tags: pl.DataFrame, reads_input: pl.DataFrame, max_errors: int
) -> pl.DataFrame:
    # First pass with exact matches
    first_pass = reads_input.join(
        parsed_tags, left_on="r2", right_on="sequence", how="left"
    )
    first_pass_mapped = first_pass.filter(~pl.col("feature_name").is_null())
    parsed_tags_extended = parsed_tags.with_columns(
        pl.col("sequence")
        .str.replace_all("A", 0)
        .str.replace_all("T", 1)
        .str.replace_all("G", 2)
        .str.replace_all("C", 3)
        .alias("ref_num")
    ).with_columns(pl.col("ref_num").str.parse_int(4))
    # Get unmapped reads
    unmapped = first_pass.filter(pl.col("feature_name").is_null())
    unmapped_with_ns = unmapped.filter(pl.col("r2").str.contains("N"))
    unmapped_without_ns = (
        unmapped.filter(~pl.col("r2").str.contains("N"))
        .with_columns(
            pl.col("r2")
            .str.replace_all("A", 0)
            .str.replace_all("T", 1)
            .str.replace_all("G", 2)
            .str.replace_all("C", 3)
            .alias("tag_num")
        )
        .with_columns(pl.col("tag_num").str.parse_int(4))
    )
    # Second pass based on closest
    second_pass = (
        unmapped_without_ns.sort("tag_num")
        .join_asof(
            parsed_tags_extended.sort("ref_num"), left_on="tag_num", right_on="ref_num"
        )
        .with_columns(diff=pl.col("ref_num") - pl.col("tag_num"))
        .select(["barcode", "umi", "r2", "feature_name_right", "diff"])
        .rename({"feature_name_right": "feature_name"})
        .with_columns(
            pl.when(pl.col("diff").abs() <= max_errors)
            .then(pl.col("feature_name"))
            .otherwise(None)
        )
    ).drop("diff")
    results = pl.concat([first_pass_mapped, second_pass, unmapped_with_ns])
    return results


def map_barcodes_pl(
    cell_reference: pl.DataFrame, mapped_reads_input: pl.DataFrame, max_errors: int
) -> pl.DataFrame:
    # First pass with exact matches
    first_pass = mapped_reads_input.join(
        cell_reference.with_row_count("barcode_id"),
        left_on="barcode",
        right_on="reference",
        how="left",
    )
    first_pass_mapped = first_pass.filter(~pl.col("barcode_id").is_null())
    ref_barcodes = first_pass_mapped.select("barcode").unique()
    ref_cells_extended = ref_barcodes.with_columns(
        pl.col("reference")
        .str.replace_all("A", 0)
        .str.replace_all("T", 1)
        .str.replace_all("G", 2)
        .str.replace_all("C", 3)
        .alias("ref_barcode_num")
    ).with_columns(pl.col("ref_barcode_num").str.parse_int(4))
    # Get unmapped reads
    unmapped = first_pass.filter(pl.col("barcode_id").is_null())
    unmapped_with_ns = unmapped.filter(pl.col("barcode").str.contains("N"))
    unmapped_without_ns = (
        unmapped.filter(~pl.col("barcode").str.contains("N"))
        .with_columns(
            pl.col("barcode")
            .str.replace_all("A", 0)
            .str.replace_all("T", 1)
            .str.replace_all("G", 2)
            .str.replace_all("C", 3)
            .alias("barcode_num")
        )
        .with_columns(pl.col("barcode_num").str.parse_int(4))
    )
    # Second pass based on closest
    second_pass = (
        unmapped_without_ns.sort("barcode_num")
        .join_asof(
            ref_cells_extended.sort("ref_barcode_num"),
            left_on="barcode_num",
            right_on="ref_barcode_num",
        )
        .with_columns(diff=pl.col("ref_barcode_num") - pl.col("barcode_num"))
        .select(["barcode", "umi", "r2", "feature_name_right", "diff"])
        .rename({"feature_name_right": "feature_name"})
        .with_columns(
            pl.when(pl.col("diff").abs() <= 1)
            .then(pl.col("feature_name"))
            .otherwise(None)
        )
    ).drop("diff")
    results = pl.concat([first_pass_mapped, second_pass, unmapped_with_ns])
    return results


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
