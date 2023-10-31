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


from cite_seq_count.preprocessing import (
    WHITELIST_COLUMN,
    BARCODE_COLUMN,
    CORRECTED_BARCODE_COLUMN,
    R2_COLUMN,
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
    barcodes_df: pl.DataFrame, barcode_subset_df: pl.DataFrame, hamming_distance: int
) -> tuple[pl.DataFrame, int]:
    print("Correcting barcodes")
    joined = (
        barcodes_df.sort(BARCODE_COLUMN)
        .join_asof(
            barcode_subset_df.sort(WHITELIST_COLUMN),
            left_on=BARCODE_COLUMN,
            right_on=WHITELIST_COLUMN,
        )
        .with_columns(
            pl.when(pl.col(BARCODE_COLUMN) == pl.col(WHITELIST_COLUMN))
            .then(False)
            .otherwise(True)
            .alias(CORRECTED_BARCODE_COLUMN)
        )
        .filter(~pl.col(WHITELIST_COLUMN).is_null())
    )
    joined = (
        pl.concat(
            [
                joined,
                joined.map_rows(
                    lambda x: distance.Hamming.distance(s1=x[0], s2=x[2])
                ).rename({"map": "hamming_distance"}),
            ],
            how="horizontal",
        )
        .filter(pl.col("hamming_distance") <= hamming_distance)
        .drop("hamming_distance")
    )
    n_corrected_barcodes = joined.filter(pl.col(CORRECTED_BARCODE_COLUMN)).shape[0]
    print("Barcodes corrected")
    return joined.drop(CORRECTED_BARCODE_COLUMN), n_corrected_barcodes


def correct_barcodes(
    barcodes_df: pl.DataFrame,
    barcode_subset: pl.DataFrame,
    collapsing_threshold: int,
):
    """Correct barcodes based on a given whitelist

    Args:
        barcodes_df (pl.DataFrame): All barcodes
        barcode_subset (pl.DataFrame): Given whitelist
        collapsing_threshold (int): Hamming distance to correct on

    Returns:
        mapped_reads_barcode_corrected (pl.DataFrame): Barcode corrected mapped reads
        n_corrected_barcodes (int): Number of corrected barcodes
    """
    print("Generating barcode tree from reference list")
    barcode_tree = pybktree.BKTree(
        Levenshtein.hamming, barcode_subset[WHITELIST_COLUMN].to_list()
    )
    print("Finding original barcodes")
    barcode_mapping = (
        barcodes_df.select(pl.col(BARCODE_COLUMN))
        .unique()
        .filter(~pl.col(BARCODE_COLUMN).is_in(barcode_subset[WHITELIST_COLUMN]))
        .with_columns(
            pl.col(BARCODE_COLUMN)
            .map_elements(
                lambda x: find_original_barcode(
                    x, barcode_tree=barcode_tree, distance=collapsing_threshold
                )
            )
            .alias(CORRECTED_BARCODE_COLUMN)
        )
        .with_columns(
            barcode_corrected=pl.col(BARCODE_COLUMN) == pl.col(CORRECTED_BARCODE_COLUMN)
        )
        .drop(BARCODE_COLUMN)
    )
    print("Collapsing wrong barcodes with original barcodes")
    mapped_reads_barcode_corrected = barcodes_df.join(
        barcode_mapping.drop(BARCODE_COLUMN),
        left_on=BARCODE_COLUMN,
        right_on=CORRECTED_BARCODE_COLUMN,
    ).drop("corrected")
    n_corrected_barcodes = barcode_mapping.filter(pl.col("corrected")).shape[0]
    return mapped_reads_barcode_corrected, n_corrected_barcodes


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
