import os

import polars as pl
from rapidfuzz import distance

from umi_tools import network

from cite_seq_count.constants import (
    SUBSET_COLUMN,
    BARCODE_COLUMN,
    FEATURE_NAME_COLUMN,
    R2_COLUMN,
    UMI_COLUMN,
    COUNT_COLUMN,
    UNMAPPED_NAME,
)


def correct_barcodes_pl(
    barcodes_df: pl.DataFrame,
    barcode_subset_df: pl.DataFrame,
    hamming_distance: int,
) -> tuple[pl.DataFrame, int, dict]:
    """Corrects barcodes using a subset based on join_asof from polars.
    Uses both forward and backward strategy to dinf the closest barcode

    Args:
        barcodes_df (pl.DataFrame): All barcodes with their respective counts
        barcode_subset_df (pl.DataFrame): Barcode reference used to correct
        hamming_distance (int): Max hamming distance allowed
        mapped_barcodes (dict): Dict of mapped barcodes

    Returns:
        tuple[pl.DataFrame, int]: The corrected version of the input barcodes_df, number of corrected barcodes
    """
    print("Correcting barcodes")
    corrected_barcodes_pl = pl.DataFrame(
        schema={
            BARCODE_COLUMN: pl.Utf8,
            "count": pl.UInt32,
            SUBSET_COLUMN: pl.Utf8,
            "hamming_distance": pl.UInt8,
        }
    )

    methods = ["forward", "backward"]
    for method in methods:
        current_barcodes = (
            barcodes_df.filter(
                (~pl.col(BARCODE_COLUMN).is_in(corrected_barcodes_pl[BARCODE_COLUMN]))
                & (~pl.col(BARCODE_COLUMN).is_in(barcode_subset_df[SUBSET_COLUMN]))
            )
            .sort(BARCODE_COLUMN)
            .join_asof(
                barcode_subset_df.sort(SUBSET_COLUMN),
                left_on=BARCODE_COLUMN,
                right_on=SUBSET_COLUMN,
                strategy=method,  # type: ignore
            )
            .filter(~pl.col(SUBSET_COLUMN).is_null())
            .with_columns(
                pl.struct(pl.col(BARCODE_COLUMN), pl.col(SUBSET_COLUMN))
                .map_elements(
                    lambda x: distance.Hamming.distance(
                        x[BARCODE_COLUMN], x[SUBSET_COLUMN]
                    ),
                    return_dtype=pl.UInt8,
                )
                .alias("hamming_distance")
            )
            .filter(pl.col("hamming_distance") <= hamming_distance)
        )
        corrected_barcodes_pl = pl.concat([corrected_barcodes_pl, current_barcodes])
    mapped_barcodes = dict(
        corrected_barcodes_pl.select(BARCODE_COLUMN, SUBSET_COLUMN).iter_rows()  # type: ignore
    )
    final_corrected = (
        barcodes_df.with_columns(
            pl.col(BARCODE_COLUMN).map_dict(mapped_barcodes, default=pl.first())
        )
        .group_by(BARCODE_COLUMN)
        .agg(pl.sum("count"))
    )
    print("Barcodes corrected")
    n_corrected_barcodes = corrected_barcodes_pl.shape[0]

    return final_corrected, n_corrected_barcodes, mapped_barcodes


def summarise_unmapped_df(main_df: pl.DataFrame, unmapped_r2_df: pl.DataFrame):
    """Merge main df and unmapped df to get a summary of the unmapped reads

    Args:
        main_df (pl.DataFrame): _description_
        unmapped_r2_df (pl.DataFrame): _description_
    """
    unmapped_r2_df = (
        unmapped_r2_df.filter(pl.col(FEATURE_NAME_COLUMN) == UNMAPPED_NAME)
        .join(main_df, on=R2_COLUMN, how="left")
        .with_columns(
            pl.when(pl.col(FEATURE_NAME_COLUMN).is_null())
            .then(pl.col(R2_COLUMN))
            .otherwise(pl.col(FEATURE_NAME_COLUMN))
            .alias(FEATURE_NAME_COLUMN)
        )
    )
    unmapped_df = unmapped_r2_df.group_by(FEATURE_NAME_COLUMN).agg(pl.count())

    return unmapped_df


def generate_umi_counts(read_counts: pl.DataFrame) -> pl.DataFrame:
    """Generate umi counts from read counts

    Args:
        read_counts (pl.DataFrame): Read counts

    Returns:
        pl.DataFrame: Umi counts
    """
    umi_counts = read_counts.group_by([BARCODE_COLUMN, FEATURE_NAME_COLUMN]).agg(
        pl.count()
    )
    return umi_counts


def update_main_df(main_df: pl.DataFrame, mapped_barcodes: dict):
    """Update the main data df with the corrected barcodes

    Args:
        main_df (pl.DataFrame): Data of all reads
        mapped_barcodes (dict): Mapped barcodes from correction

    Returns:
        pl.DataFrame: Data of all reads with barcodes corrected
    """
    main_df = (
        main_df.with_columns(
            pl.col(BARCODE_COLUMN).map_dict(mapped_barcodes, default=pl.first())
        )
        .group_by([BARCODE_COLUMN, UMI_COLUMN, R2_COLUMN])
        .agg(pl.sum(COUNT_COLUMN))
    )
    return main_df


def correct_umis_df(
    main_df, mapped_r2_df, umi_distance=1, cluster_method="directional", max_umis=20000
):
    merged = mapped_r2_df.join(main_df, on=R2_COLUMN)
    temp = (
        merged.with_columns(umi=pl.col(UMI_COLUMN).cast(pl.Binary))
        .group_by([R2_COLUMN, FEATURE_NAME_COLUMN, BARCODE_COLUMN])
        .agg(pl.struct(pl.col(UMI_COLUMN), pl.col(COUNT_COLUMN)))
    )
    clustered_cells = (
        temp.filter(pl.col(UMI_COLUMN).list.len() > max_umis)
        .select(BARCODE_COLUMN)
        .unique()
        .get_column(BARCODE_COLUMN)
        .to_list()
    )

    mapping_list = []
    umi_clusterer = network.UMIClusterer(cluster_method=cluster_method)
    for r2, feature_name, barcode, umis in temp.filter(
        (pl.col(UMI_COLUMN).list.len() > 1) & (pl.col(UMI_COLUMN).list.len() < max_umis)
    ).iter_rows():
        corrected_umis = correct_umis(
            umis,
            umi_distance=umi_distance,
            umi_clusterer=umi_clusterer,
        )
        if len(corrected_umis) == 0:
            continue
        for umi_set in corrected_umis:
            for index, umi in enumerate(umi_set):
                if index != 0:
                    mapping_list.append([r2, feature_name, barcode, umi, umi_set[0]])
    mapping_df = (
        pl.DataFrame(
            mapping_list,
            schema={
                R2_COLUMN: pl.String,
                FEATURE_NAME_COLUMN: pl.String,
                BARCODE_COLUMN: pl.String,
                "orig": pl.Binary,
                "replace": pl.Binary,
            },
        )
        .with_columns(
            pl.col("orig").cast(pl.String).alias(UMI_COLUMN),
            pl.col("replace").cast(pl.String),
        )
        .drop("orig")
    )
    read_counts = (
        merged.join(
            mapping_df,
            on=[R2_COLUMN, FEATURE_NAME_COLUMN, BARCODE_COLUMN, UMI_COLUMN],
            how="left",
        )
        .with_columns(
            pl.when(pl.col("replace").is_null())
            .then(pl.col(UMI_COLUMN))
            .otherwise(pl.col("replace"))
            .alias(UMI_COLUMN)
        )
        .drop("replace")
        .group_by([R2_COLUMN, BARCODE_COLUMN, FEATURE_NAME_COLUMN, UMI_COLUMN])
        .agg(pl.sum(COUNT_COLUMN))
        .drop(R2_COLUMN)
    )
    n_corrected_umis = mapping_df.shape[0]
    return read_counts, n_corrected_umis, clustered_cells


def correct_umis(umis_list, umi_distance, umi_clusterer):
    umis = dict([(i["umi"], i["count"]) for i in umis_list])

    res = umi_clusterer(umis, umi_distance)
    corrected = [corrected_umis for corrected_umis in res if len(corrected_umis) > 1]
    return corrected


# def generate_mtx_counts(
#     main_df: pl.DataFrame,
#     barcode_subset: pl.DataFrame,
#     mapped_r2_df: pl.DataFrame,
#     data_type: str,
# ) -> pl.DataFrame:
#     if data_type == "read":
#         return (
#             main_df.join(barcode_subset, on=BARCODE_COLUMN)
#             .join(mapped_r2_df, on=R2_COLUMN)
#             .group_by([BARCODE_COLUMN, FEATURE_NAME_COLUMN])
#             .agg(pl.sum(COUNT_COLUMN))
#         )
#     else:
#         return (
#             main_df.join(barcode_subset, on=BARCODE_COLUMN)
#             .join(mapped_r2_df, on=R2_COLUMN)
#             .group_by([BARCODE_COLUMN, FEATURE_NAME_COLUMN])
#             .agg(pl.count())
#         )
