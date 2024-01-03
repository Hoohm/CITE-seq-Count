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
        .agg(pl.sum("count"))
    )
    return main_df


def find_corrected_umis(
    main_df, mapped_r2_df, umi_distance=1, cluster_method="directional", max_umis=20000
):
    merged = mapped_r2_df.join(main_df, on="r2")
    temp = (
        merged.with_columns(umi=pl.col("umi").cast(pl.Binary))
        .group_by(["r2", "feature_name", "barcode"])
        .agg(pl.struct(pl.col("umi"), pl.col("count")))
        .filter((pl.col("umi").list.len() > 1) & (pl.col("umi").list.len() < max_umis))
    )
    mapping = {"r2": [], "feature_name": [], "barcode": [], "orig": [], "replace": []}
    for r2, feature_name, barcode, umis in temp.iter_rows():
        corrected_umis = correct_umis(
            umis, cluster_method=cluster_method, umi_distance=umi_distance
        )
        if len(corrected_umis) == 0:
            continue
        for umi_set in corrected_umis:
            for index, umi in enumerate(umi_set):
                if index != 0:
                    mapping["r2"].append(r2)
                    mapping["feature_name"].append(feature_name)
                    mapping["barcode"].append(barcode)
                    mapping["orig"].append(umi)
                    mapping["replace"].append(umi_set[0])
    mapping_df = (
        pl.DataFrame(mapping)
        .with_columns(
            pl.col("orig").cast(pl.String).alias("umi"),
            pl.col("replace").cast(pl.String),
        )
        .drop("orig")
    )
    read_counts = (
        merged.join(mapping_df, on=["r2", "feature_name", "barcode", "umi"], how="left")
        .with_columns(
            pl.when(pl.col("replace").is_null())
            .then(pl.col("umi"))
            .otherwise(pl.col("replace"))
            .alias("umi")
        )
        .drop("replace")
        .group_by(["r2", "barcode", "feature_name", "umi"])
        .agg(pl.sum("count")).drop("r2")
    )
    return read_counts


def correct_umis(umis_list, cluster_method, umi_distance):
    umi_clusterer = network.UMIClusterer(cluster_method=cluster_method)
    umis = dict([(i["umi"], i["count"]) for i in umis_list])

    res = umi_clusterer(umis, umi_distance)
    corrected = [corrected_umis for corrected_umis in res if len(corrected_umis) > 1]
    return corrected


# UMI correction section
def correct_umis_in_cells(umi_correction_input):
    """
    Corrects umi barcodes within same cell/tag groups.

    Args:
        final_results (dict): Dict of dict of Counters with mapping results.
        collapsing_threshold (int): Max distance between umis.
        filtered_cells (set): Set of cells to go through.
        max_umis (int): Maximum UMIs to consider for one cluster.

    Returns:
        final_results (dict): Same as input but with corrected umis.
        corrected_umis (int): How many umis have been corrected.
        clustered_umi_count_cells (set): Set of uncorrected cells.
    """

    (final_results, collapsing_threshold, max_umis, unmapped_id) = umi_correction_input
    print(
        "Started umi correction in child process {} working on {} cells".format(
            os.getpid(), len(final_results)
        )
    )
    corrected_umis = 0
    clustered_cells = set()
    cells = final_results.keys()
    for cell_barcode in cells:
        for TAG in final_results[cell_barcode]:
            if TAG == unmapped_id:
                final_results[cell_barcode].pop(unmapped_id)

            n_umis = len(final_results[cell_barcode][TAG])
            if n_umis > 1 and n_umis <= max_umis:
                umi_clusters = network.UMIClusterer()
                UMIclusters = umi_clusters(
                    final_results[cell_barcode][TAG], collapsing_threshold
                )
                (new_res, temp_corrected_umis) = update_umi_counts(
                    UMIclusters, final_results[cell_barcode][TAG]
                )
                final_results[cell_barcode][TAG] = new_res
                corrected_umis += temp_corrected_umis
            elif n_umis > max_umis:
                clustered_cells.add(cell_barcode)
    print(f"Finished correcting umis in child {os.getpid()}")
    return (final_results, corrected_umis, clustered_cells)


def generate_mtx_counts(
    main_df: pl.DataFrame,
    barcode_subset: pl.DataFrame,
    mapped_r2_df: pl.DataFrame,
    data_type: str,
) -> pl.DataFrame:
    if data_type == "read":
        return (
            main_df.join(barcode_subset, on=BARCODE_COLUMN)
            .join(mapped_r2_df, on=R2_COLUMN)
            .group_by([BARCODE_COLUMN, FEATURE_NAME_COLUMN])
            .agg(pl.sum(COUNT_COLUMN))
        )
    else:
        return (
            main_df.join(barcode_subset, on=BARCODE_COLUMN)
            .join(mapped_r2_df, on=R2_COLUMN)
            .group_by([BARCODE_COLUMN, FEATURE_NAME_COLUMN])
            .agg(pl.count())
        )
