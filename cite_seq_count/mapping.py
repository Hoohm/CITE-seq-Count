import time
import csv

import sys
import os
import Levenshtein

from collections import Counter
from collections import defaultdict
from collections import namedtuple

# pylint: disable=no-name-in-module
from multiprocess import Pool

from cite_seq_count.processing import merge_results
from cite_seq_count import secondsToText


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
    """Read through R1/R2 files and generate.

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
