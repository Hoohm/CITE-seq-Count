import os
import Levenshtein
import pybktree

import polars as pl
from rapidfuzz import distance

from collections import namedtuple

# pylint: disable=no-name-in-module
from multiprocess import Pool


from numpy import int32
from scipy import sparse
from umi_tools import network


from cite_seq_count.constants import (
    WHITELIST_COLUMN,
    BARCODE_COLUMN,
    CORRECTED_BARCODE_COLUMN,
    R2_COLUMN,
    UMI_COLUMN,
)

# Unit Barcode correction


def find_original_barcode(barcode: str, barcode_tree: pybktree.BKTree, distance: int):
    """Pare a BKtree to find the original barcode to correct to.

    Args:
        barcode (str): barcode to be corrected
        barcode_tree (pybktree.BKTree): Barcode whitelist BKTree
        distance (int): Hamming distance to search for

    Returns:
        barcode(str): corrected barcode
    """
    candidates = [
        white_cell for d, white_cell in barcode_tree.find(barcode, distance) if d > 0
    ]
    if len(candidates) == 1:
        barcode = candidates[0]
    return barcode


def merge_results(
    mapped_r2_df: pl.DataFrame,
    corrected_barcodes_df: pl.DataFrame,
    input_df: pl.DataFrame,
):
    merged = (
        input_df.join(mapped_r2_df, on=R2_COLUMN, how="inner")
        .join(corrected_barcodes_df.drop("count"), on=BARCODE_COLUMN, how="inner")
        .drop([R2_COLUMN, BARCODE_COLUMN])
    )
    return merged


def correct_barcodes_pl(
    barcodes_df: pl.DataFrame,
    barcode_subset_df: pl.DataFrame,
    hamming_distance: int,
) -> tuple[pl.DataFrame, int, dict]:
    """Corrects barcodes using a whitelist based on join_asof from polars.
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
            "count": pl.Int64,
            WHITELIST_COLUMN: pl.Utf8,
            "hamming_distance": pl.UInt8,
        }
    )
    methods = ["backward", "forward"]
    for method in methods:
        current_barcodes = (
            barcodes_df.filter(
                (~pl.col(BARCODE_COLUMN).is_in(corrected_barcodes_pl[BARCODE_COLUMN]))
                & (~pl.col(BARCODE_COLUMN).is_in(barcode_subset_df[WHITELIST_COLUMN]))
            )
            .sort(BARCODE_COLUMN)
            .join_asof(
                barcode_subset_df.sort(WHITELIST_COLUMN),
                left_on=BARCODE_COLUMN,
                right_on=WHITELIST_COLUMN,
                strategy=method,
            )
            .filter(~pl.col(WHITELIST_COLUMN).is_null())
            .with_columns(
                pl.struct(pl.col(BARCODE_COLUMN), pl.col(WHITELIST_COLUMN))
                .map_elements(
                    lambda x: distance.Hamming.distance(
                        x[BARCODE_COLUMN], x[WHITELIST_COLUMN]
                    ),
                    return_dtype=pl.UInt8,
                )
                .alias("hamming_distance")
            )
            .filter(pl.col("hamming_distance") <= hamming_distance)
        )
        corrected_barcodes_pl = pl.concat([corrected_barcodes_pl, current_barcodes])
    mapped_barcodes = dict(
        corrected_barcodes_pl.select(BARCODE_COLUMN, WHITELIST_COLUMN).iter_rows()
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


def update_main_df(main_df: pl.DataFrame, mapped_barcodes: dict):
    """Update the main data df with the corrected barcodes

    Args:
        main_df (pl.DataFrame): Data of all reads
        mapped_barcodes (dict): Mapped barcodes from correction

    Returns:
        pl.DataFrame: Data of all reads with barcodes corrected
    """
    main_df = (
        main_df.with_columns(pl.col(BARCODE_COLUMN).map_dict(mapped_barcodes))
        .group_by([BARCODE_COLUMN, UMI_COLUMN, R2_COLUMN])
        .agg(pl.sum("count"))
    )
    return main_df


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


def update_umi_counts(UMIclusters, cell_tag_counts):
    """
    Update a dict object with umis corrected.

    Args:
        UMIclusters (list): List of lists with corrected umis
        cell_tag_counts (Counter): Counter of umis

    Returns:
        cell_tag_counts (Counter): Updated Counter of umis
        temp_corrected_umis (int): Number of corrected umis
    """
    temp_corrected_umis = 0
    for (
        umi_cluster
    ) in UMIclusters:  # This is a list with the first element the dominant barcode
        if len(umi_cluster) > 1:  # This means we got a correction
            major_umi = umi_cluster[0]
            for minor_umi in umi_cluster[1:]:
                temp_corrected_umis += 1
                temp = cell_tag_counts.pop(minor_umi)
                cell_tag_counts[major_umi] += temp
    return (cell_tag_counts, temp_corrected_umis)


def run_umi_correction(final_results, filtered_cells, unmapped_id, args):
    input_queue = []
    umi_correction_input = namedtuple(
        "umi_correction_input",
        ["cells", "collapsing_threshold", "max_umis", "unmapped_id"],
    )
    cells_results = {}
    n_cells = 0
    num_chunks = 0

    print("preparing UMI correction jobs")
    cell_batch_size = round(len(filtered_cells) / args.n_threads) + 1
    for cell in filtered_cells:
        cells_results[cell] = final_results.pop(cell)
        n_cells += 1
        if n_cells % cell_batch_size == 0:
            input_queue.append(
                umi_correction_input(
                    cells=cells_results,
                    collapsing_threshold=args.umi_threshold,
                    max_umis=20000,
                    unmapped_id=unmapped_id,
                )
            )
            cells_results = {}
            num_chunks += 1

    del final_results

    input_queue.append(
        umi_correction_input(
            cells=cells_results,
            collapsing_threshold=args.umi_threshold,
            max_umis=20000,
            unmapped_id=unmapped_id,
        )
    )
    parallel_results = []
    if args.n_threads != 1:
        pool = Pool(processes=args.n_threads)
        errors = []
        correct_umis = pool.map_async(
            correct_umis_in_cells,
            input_queue,
            callback=parallel_results.append,
            error_callback=errors.append,
        )

        correct_umis.wait()
        pool.close()
        pool.join()

        if len(errors) != 0:
            for error in errors:
                print("There was an error {}", error)
    else:
        single_thread_result = correct_umis_in_cells(input_queue[0])
        parallel_results.append([single_thread_result])
    final_results = {}
    umis_corrected = 0
    clustered_cells = set()
    for chunk in parallel_results[0]:
        (temp_results, temp_umis, temp_clustered_cells) = chunk
        final_results.update(temp_results)
        umis_corrected += temp_umis
        clustered_cells.update(temp_clustered_cells)

    return final_results, umis_corrected, clustered_cells


def generate_sparse_matrices(
    final_results, parsed_tags, filtered_cells, umi_counts=False
):
    """
    Create two sparse matrices with umi and read counts.

    Args:
        final_results (dict): Results in a dict of dicts of Counters.
        parsed_tags (list): Ordered tags in a list of tuples.

    Returns:
        results_matrix (scipy.sparse.dok_matrix): UMI or Read counts


    """
    unmapped_id = len(parsed_tags)
    if umi_counts:
        n_features = len(parsed_tags)
    else:
        n_features = len(parsed_tags) + 1
    results_matrix = sparse.dok_matrix((n_features, len(filtered_cells)), dtype=int32)

    for i, cell_barcode in enumerate(filtered_cells):
        if cell_barcode not in final_results.keys():
            continue
        for TAG_id in final_results[cell_barcode]:
            # if TAG_id in final_results[cell_barcode]:
            if umi_counts:
                if TAG_id == unmapped_id:
                    continue
                else:
                    results_matrix[TAG_id, i] = len(final_results[cell_barcode][TAG_id])
            else:
                results_matrix[TAG_id, i] = sum(
                    final_results[cell_barcode][TAG_id].values()
                )
    return results_matrix
