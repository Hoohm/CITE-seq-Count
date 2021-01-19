import time
import sys
import os
import Levenshtein
import pybktree
import csv

from collections import Counter
from collections import defaultdict
from collections import namedtuple

# pylint: disable=no-name-in-module
from multiprocess import Pool

from itertools import islice
from numpy import int32
from scipy import sparse
from umi_tools import network
import umi_tools.whitelist_methods as whitelist_methods


from cite_seq_count import secondsToText


def find_best_match(TAG_seq, tags, maximum_distance):
    """
    Find the best match from the list of tags.

    Compares the Levenshtein distance between tags and the trimmed sequences.
    The tag and the sequence must have the same length.
    If no matches found returns 'unmapped'.
    We add 1
    Args:
        TAG_seq (string): Sequence from R2 already start trimmed
        tags (dict): A dictionary with the TAGs as keys and TAG Names as values.
        maximum_distance (int): Maximum distance given by the user.

    Returns:
        best_match (string): The TAG name that will be used for counting.
    """
    best_match = len(tags)
    best_score = maximum_distance
    for tag in tags:
        # pylint: disable=no-member
        score = Levenshtein.hamming(tag.sequence, TAG_seq[: len(tag.sequence)])
        if score == 0:
            # Best possible match
            return tag.id
        elif score <= best_score:
            best_score = score
            best_match = tag.id
            return best_match
    return best_match


def find_best_match_shift(TAG_seq, tags):
    """
    Find the best match from the list of tags with sliding window.
    Only works with exact match.
    Just checks if the string is in the sequence.
    If no matches found returns 'unmapped'.
    
    Args:
        TAG_seq (string): Sequence from R2 already start trimmed
        tags (dict): A dictionary with the TAGs as keys and TAG Names as values.

    Returns:
        best_match (string): The TAG name that will be used for counting.
    """
    best_match = "unmapped"
    for tag in tags:
        if tag.sequence in TAG_seq:
            return tag.name
    return best_match


def map_reads(mapping_input):
    """Read through R1/R2 files and generate a islice starting at a specific index.

    It reads both Read1 and Read2 files, creating a dict based on cell barcode.

    Args:
        mapping_input (namedtuple): List of paramters to run in parallel.
            filename (str): Path to the chunk file
            tags (list): List of named tuples tags
            debug (bool): Should debug information be shown or not
            maximum_distance (int): Maximum distance given by the user
            sliding_window (bool): A bool enabling a sliding window search

    Returns:
        results (dict): A dict of dict of Counters with the mapping results.
        no_match (Counter): A counter with unmapped sequences.
    """
    # Initiate values
    (filename, tags, debug, maximum_distance, sliding_window) = mapping_input
    print("Started mapping in child process {}".format(os.getpid()))
    results = {}
    no_match = Counter()
    n = 1

    unmapped_id = len(tags)
    # Progress info
    t = time.time()
    with open(filename, "r") as input_file:
        reads = csv.reader(input_file)
        for read in reads:
            cell_barcode = read[0]
            # This change in bytes is required by umi_tools for umi correction
            UMI = bytes(read[1], "ascii")
            read2 = read[2]
            if n % 1000000 == 0:
                print(
                    "Processed 1,000,000 reads in {}. Total "
                    "reads: {:,} in child {}".format(
                        secondsToText.secondsToText(time.time() - t), n, os.getpid()
                    )
                )
                sys.stdout.flush()
                t = time.time()

            if cell_barcode not in results:
                results[cell_barcode] = defaultdict(Counter)

            if sliding_window:
                best_match = find_best_match_shift(read2, tags)
            else:
                best_match = find_best_match(read2, tags, maximum_distance)

            results[cell_barcode][best_match][UMI] += 1

            if best_match == unmapped_id:
                no_match[read2] += 1

            if debug:
                print(
                    "cell_barcode:{0}\tUMI:{1}\tTAG_seq:{2}\n"
                    "cell barcode length:{3}\tUMI length:{4}\tTAG sequence length:{5}\n"
                    "Best match is: {6}\n".format(
                        cell_barcode,
                        UMI,
                        read2,
                        len(cell_barcode),
                        len(UMI),
                        len(read2),
                        tags[best_match].name,
                    )
                )
                sys.stdout.flush()
            n += 1
        print(
            "Mapping done for process {}. Processed {:,} reads".format(
                os.getpid(), n - 1
            )
        )
        sys.stdout.flush()

    return (results, no_match)


def merge_results(parallel_results, unmapped_id):
    """Merge chunked results from parallel processing.

    Args:
        parallel_results (list): List of dict with mapping results.

    Returns:
        merged_results (dict): Results combined as a dict of dict of Counters
        umis_per_cell (Counter): Total umis per cell as a Counter
        reads_per_cell (Counter): Total reads per cell as a Counter
        merged_no_match (Counter): Unmapped tags as a Counter
    """
    merged_results = {}
    merged_no_match = Counter()
    umis_per_cell = Counter()
    reads_per_cell = Counter()
    for chunk in parallel_results:
        mapped = chunk[0]
        unmapped = chunk[1]
        for cell_barcode in mapped:
            if cell_barcode not in merged_results:
                merged_results[cell_barcode] = defaultdict(Counter)
            for TAG in mapped[cell_barcode]:
                # We don't want to capture unmapped data in the umi counts
                if TAG == unmapped_id:
                    continue
                # Test the counter. Returns false if empty
                if mapped[cell_barcode][TAG]:
                    for UMI in mapped[cell_barcode][TAG]:
                        merged_results[cell_barcode][TAG][UMI] += mapped[cell_barcode][
                            TAG
                        ][UMI]
                        umis_per_cell[cell_barcode] += len(mapped[cell_barcode][TAG])
                        reads_per_cell[cell_barcode] += mapped[cell_barcode][TAG][UMI]
        merged_no_match.update(unmapped)
    return merged_results, umis_per_cell, reads_per_cell, merged_no_match


def check_unmapped(no_match, too_short, total_reads, start_trim):
    """Check if the number of unmapped is higher than 99%"""
    sum_unmapped = sum(no_match.values()) + too_short
    if sum_unmapped / total_reads > float(0.99):
        sys.exit(
            """More than 99% of your data is unmapped.\nPlease check that your --start_trim {} parameter is correct and that your tags file is properly formatted""".format(
                start_trim
            )
        )


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
    print("Finished correcting umis in child {}".format(os.getpid()))
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


def collapse_cells(true_to_false, umis_per_cell, final_results, ab_map):
    """
    Collapses cell barcodes based on the mapping true_to_false

    Args:
        true_to_false (dict): Mapping between the reference and the "mutated" barcodes.
        umis_per_cell (Counter): Counter of number of umis per cell.
        final_results (dict): Dict of dict of Counters with mapping results.
        ab_map (dict): Dict of the TAGS.

    Returns:
        umis_per_cell (Counter): Counter of number of umis per cell.
        final_results (dict): Same as input but with corrected cell barcodes.
        corrected_barcodes (int): How many cell barcodes have been corrected.
    """
    print("Collapsing cell barcodes")
    corrected_barcodes = 0
    for real_barcode in true_to_false:
        # If the cell barcode is not in the results
        # add it in.
        if real_barcode not in final_results:
            final_results[real_barcode] = defaultdict()
            for TAG in ab_map:
                final_results[real_barcode][TAG.id] = Counter()
        for wrong_barcode in true_to_false[real_barcode]:
            temp = final_results.pop(wrong_barcode)

            for TAG in temp.keys():
                if TAG in final_results[real_barcode]:
                    final_results[real_barcode][TAG].update(temp[TAG])
                else:
                    final_results[real_barcode][TAG] = temp[TAG]
            corrected_barcodes += 1
            temp_umi_counts = umis_per_cell.pop(wrong_barcode)
            # temp_read_counts = reads_per_cell.pop(wrong_barcode)

            umis_per_cell[real_barcode] += temp_umi_counts
            # reads_per_cell[real_barcode] += temp_read_counts

    return (umis_per_cell, final_results, corrected_barcodes)


def correct_cells_no_reference_list(
    final_results,
    reads_per_cell,
    umis_per_cell,
    collapsing_threshold,
    expected_cells,
    ab_map,
):
    """
    Corrects cell barcodes.
    
    Args:
        final_results (dict): Dict of dict of Counters with mapping results.
        umis_per_cell (Counter): Counter of number of umis per cell.
        collapsing_threshold (int): Max distance between umis.
        expected_cells (int): Number of expected cells.
        ab_map (dict): Dict of the TAGS.
    
    Returns:
        final_results (dict): Same as input but with corrected umis.
        umis_per_cell (Counter): Counter of umis per cell after cell barcode correction
        corrected_umis (int): How many umis have been corrected.
    """
    print("Looking for a reference list")
    _, true_to_false = whitelist_methods.getCellWhitelist(
        knee_method="density",
        cell_barcode_counts=reads_per_cell,
        expect_cells=expected_cells,
        cell_number=expected_cells,
        error_correct_threshold=collapsing_threshold,
        plotfile_prefix=False,
    )
    if true_to_false is None:
        print("Failed to find a good reference list. Will not correct cell barcodes")
        corrected_barcodes = 0
        return (final_results, umis_per_cell, corrected_barcodes)
    (umis_per_cell, final_results, corrected_barcodes) = collapse_cells(
        true_to_false=true_to_false,
        umis_per_cell=umis_per_cell,
        final_results=final_results,
        ab_map=ab_map,
    )
    return (final_results, umis_per_cell, corrected_barcodes)


def correct_cells_reference_list(
    final_results, umis_per_cell, reference_list, collapsing_threshold, ab_map
):
    """
    Corrects cell barcodes based on a given reference_list.
    
    Args:
        final_results (dict): Dict of dict of Counters with mapping results.
        umis_per_cell (Counter): Counter of UMIs per cell.
        reference_list (set): The reference_list reference given by the user.
        collapsing_threshold (int): Max distance between umis.
        ab_map (OrederedDict): Tags in an ordered dict.

    
    Returns:
        final_results (dict): Same as input but with corrected umis.
        umis_per_cell (Counter): Updated UMI counts after correction.
        corrected_barcodes (int): How many umis have been corrected.
    """
    print("Generating barcode tree from reference list")
    # pylint: disable=no-member
    barcode_tree = pybktree.BKTree(Levenshtein.hamming, reference_list)
    barcodes = set(umis_per_cell)
    print("Selecting reference candidates")
    print("Processing {:,} cell barcodes".format(len(barcodes)))

    # Run with one process
    true_to_false = find_true_to_false_map(
        barcode_tree=barcode_tree,
        cell_barcodes=barcodes,
        reference_list=reference_list,
        collapsing_threshold=collapsing_threshold,
    )
    print("Collapsing wrong barcodes with original barcodes")
    (umis_per_cell, final_results, corrected_barcodes) = collapse_cells(
        true_to_false, umis_per_cell, final_results, ab_map
    )
    return (final_results, umis_per_cell, corrected_barcodes)


def find_true_to_false_map(
    barcode_tree, cell_barcodes, reference_list, collapsing_threshold
):
    """
    Creates a mapping between "fake" cell barcodes and their original true barcode.

    Args:
        barcode_tree (BKTree): BKTree of all original cell barcodes.
        cell_barcodes (List): Cell barcodes to go through.
        reference_list (dict): Dict of the reference_list, the "true" cell barcodes.
        collasping_threshold (int): How many mistakes to correct.

    Return:
        true_to_false (defaultdict(list)): Contains the mapping between the fake and real barcodes. The key is the real one.
    """
    true_to_false = defaultdict(list)
    for cell_barcode in cell_barcodes:
        if cell_barcode in reference_list:
            # if the barcode is already reference_listed, no need to add
            continue
        # get all members of reference_list that are at distance of collapsing_threshold
        candidates = [
            white_cell
            for d, white_cell in barcode_tree.find(cell_barcode, collapsing_threshold)
            if d > 0
        ]
        if len(candidates) == 1:
            white_cell_str = candidates[0]
            true_to_false[white_cell_str].append(cell_barcode)
        else:
            # the cell doesnt match to any reference_listed barcode,
            # hence we have to drop it
            # (as it cannot be asscociated with any frequent barcode)
            continue
    return true_to_false


def generate_sparse_matrices(
    final_results, ordered_tags, filtered_cells, umi_counts=False
):
    """
    Create two sparse matrices with umi and read counts.

    Args:
        final_results (dict): Results in a dict of dicts of Counters.
        ordered_tags (list): Ordered tags in a list of tuples.

    Returns:
        results_matrix (scipy.sparse.dok_matrix): UMI counts


    """
    unmapped_id = len(ordered_tags)
    if umi_counts:
        n_features = len(ordered_tags)
    else:
        n_features = len(ordered_tags) + 1
    results_matrix = sparse.dok_matrix((n_features, len(filtered_cells)), dtype=int32)
    # print(ordered_tags)

    for i, cell_barcode in enumerate(filtered_cells):
        if cell_barcode not in final_results.keys():
            continue
        for TAG_id in final_results[cell_barcode]:
            # if TAG_id in final_results[cell_barcode]:
            if umi_counts:
                if TAG_id == unmapped_id:
                    continue
                results_matrix[TAG_id, i] = len(final_results[cell_barcode][TAG_id])
            else:
                results_matrix[TAG_id, i] = sum(
                    final_results[cell_barcode][TAG_id].values()
                )
    return results_matrix


def map_data(input_queue, unmapped_id, args):
    """
    Maps the data given an input_queue

    Args:
        input_queue (list): List of parameters to run in parallel
        args (argparse): List of arguments
    
    Returns:
        final_results (dict): final dictionnary with results
        umis_per_cell (Counter): Counter of UMIs per cell
        reads_per_cell (Counter): Counter of reads per cell
        merged_no_match (Counter): Counter of unmapped reads
    """
    # Initialize the counts dicts that will be generated from each input fastq pair
    final_results = defaultdict(lambda: defaultdict(Counter))
    umis_per_cell = Counter()
    reads_per_cell = Counter()
    merged_no_match = Counter()

    print("Started mapping")
    parallel_results = []

    if args.n_threads == 1:
        mapped_reads = map_reads(input_queue[0])
        parallel_results.append([mapped_reads])
    else:
        pool = Pool(processes=args.n_threads)
        errors = []
        mapping = pool.map_async(
            map_reads,
            input_queue,
            callback=parallel_results.append,
            error_callback=errors.append,
        )
        mapping.wait()

        pool.close()
        pool.join()
        if len(errors) != 0:
            for error in errors:
                print(error)

        print("Merging results")
    (final_results, umis_per_cell, reads_per_cell, merged_no_match,) = merge_results(
        parallel_results=parallel_results[0], unmapped_id=unmapped_id
    )

    return final_results, umis_per_cell, reads_per_cell, merged_no_match


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


def run_cell_barcode_correction(
    final_results, umis_per_cell, reads_per_cell, reference_dict, ordered_tags, args
):
    if len(umis_per_cell) <= args.expected_cells:
        print(
            "Number of expected cells, {}, is higher "
            "than number of cells found {}.\nNot performing "
            "cell barcode correction"
            "".format(args.expected_cells, len(umis_per_cell))
        )
        bcs_corrected = 0
    else:
        print("Correcting cell barcodes")
        if not reference_dict:
            print("Reference list not given")
            (
                final_results,
                umis_per_cell,
                bcs_corrected,
            ) = correct_cells_no_reference_list(
                final_results=final_results,
                reads_per_cell=reads_per_cell,
                umis_per_cell=umis_per_cell,
                expected_cells=args.expected_cells,
                collapsing_threshold=args.bc_threshold,
                ab_map=ordered_tags,
            )
        else:
            print("Reference list given")
            (
                final_results,
                umis_per_cell,
                bcs_corrected,
            ) = correct_cells_reference_list(
                final_results=final_results,
                umis_per_cell=umis_per_cell,
                reference_list=set(reference_dict.keys()),
                collapsing_threshold=args.bc_threshold,
                ab_map=ordered_tags,
            )
    return final_results, umis_per_cell, bcs_corrected
