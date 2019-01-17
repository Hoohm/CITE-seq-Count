import time
import gzip
import sys
import os
import Levenshtein
import regex

from collections import Counter
from collections import defaultdict

from itertools import islice
from numpy import int16
from scipy import sparse
from umi_tools import network
from umi_tools import umi_methods


from cite_seq_count import secondsToText


def find_best_match(TAG_seq, tags, maximum_distance):
    """
    Find the best match from the list of tags.

    Compares the Levenshtein distance between tags and the trimmed sequences.
    The tag and the sequence must have the same length.
    If no matches found returns 'unmapped'.
    We add 1
    Args:
        TAG_seq (string): Sequence from R1 already start trimmed
        tags (dict): A dictionary with the TAGs as keys and TAG Names as values.
        maximum_distance (int): Maximum distance given by the user.

    Returns:
        best_match (string): The TAG name that will be used for counting.
    """
    best_match = 'unmapped'
    best_score = maximum_distance
    for tag, name in tags.items():
        score = Levenshtein.distance(tag, TAG_seq[:len(tag)])
        if score <= best_score:
            best_score = score
            best_match = name
            return(best_match)
    return(best_match)


def map_reads(read1_path, read2_path, tags, barcode_slice,
                umi_slice, indexes, whitelist, debug,
                start_trim, maximum_distance):
    """Read through R1/R2 files and generate a islice starting at a specific index.

    It reads both Read1 and Read2 files, creating a dict based on cell barcode.

    Args:
        read1_path (string): Path to R1.fastq.gz
        read2_path (string): Path to R2.fastq.gz
        chunk_size (int): The number of lines to process 
        tags (dict): A dictionary with the TAGs + TAG Names.
        barcode_slice (slice): A slice for extracting the Barcode portion from the
            sequence.
        umi_slice (slice): A slice for extracting the UMI portion from the
            sequence.
        indexes (list): Pair of first and last index for islice
        whitelist (set): The set of white-listed barcodes.
        debug (bool): Print debug messages. Default is False.
        start_trim (int): Number of bases to trim at the start.
        maximum_distance (int): Maximum distance given by the user.

    Returns:
        results (dict): A dict of dict of Counters with the mapping results.
        no_match (Counter): A counter with unmapped sequences.
    """
    # Initiate values
    results = {}
    no_match = Counter()
    n = 1
    t = time.time()
    with gzip.open(read1_path, 'rt') as textfile1, \
         gzip.open(read2_path, 'rt') as textfile2:
        
        # Read all 2nd lines from 4 line chunks. If first_n not None read only 4 times the given amount.
        secondlines = islice(zip(textfile1, textfile2), indexes[0]*4+1, indexes[1]*4+1, 4)
        for read1, read2 in secondlines:
            read1 = read1.strip()
            read2 = read2.strip()
            
            # Progress info
            if n % 1000000 == 0:
                print("Processed 1,000,000 reads in {}. Total "
                    "reads: {:,} in child {}".format(
                        secondsToText.secondsToText(time.time() - t),
                        n,
                        os.getpid())
                    )
                sys.stdout.flush()
                t = time.time()

            # Get cell and umi barcodes.
            cell_barcode = read1[barcode_slice]
            # This change in bytes is required by umi_tools for umi correction
            UMI = bytes(read1[umi_slice], 'ascii')
            # Trim potential starting sequences
            TAG_seq = read2[start_trim:]

            if cell_barcode not in results:
                results[cell_barcode] = defaultdict(Counter)
            
            best_match = find_best_match(TAG_seq, tags, maximum_distance)
            
            results[cell_barcode][best_match][UMI] += 1
            
            if(best_match == 'unmapped'):
                no_match[TAG_seq] += 1 
            
            if debug:
                print(
                    "\nline:{0}\n"
                    "cell_barcode:{1}\tUMI:{2}\tTAG_seq:{3}\n"
                    "line length:{4}\tcell barcode length:{5}\tUMI length:{6}\tTAG sequence length:{7}\n"
                    "Best match is: {8}"
                    .format(read1 + read2, cell_barcode, UMI, TAG_seq,
                            len(read1 + read2), len(cell_barcode), len(UMI), len(TAG_seq), best_match
                    )
                )
                sys.stdout.flush()
            n += 1
    print("Mapping done for process {}. Processed {:,} reads".format(os.getpid(), n - 1))
    sys.stdout.flush()
    return(results, no_match)


def merge_results(parallel_results):
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
                # Test the counter. Returns false if empty
                if mapped[cell_barcode][TAG]:
                    for UMI in mapped[cell_barcode][TAG]:
                        merged_results[cell_barcode][TAG][UMI] += mapped[cell_barcode][TAG][UMI]
                        if TAG != 'unmapped':
                            umis_per_cell[cell_barcode] += len(mapped[cell_barcode][TAG])
                            reads_per_cell[cell_barcode] += mapped[cell_barcode][TAG][UMI]
        merged_no_match.update(unmapped)
    return(merged_results, umis_per_cell, reads_per_cell, merged_no_match)


def correct_umis(final_results, collapsing_threshold):
    """
    Corrects umi barcodes within same cell/tag groups.
    
    Args:
        final_results (dict): Dict of dict of Counters with mapping results.
        collapsing_threshold (int): Max distance between umis.
    
    Returns:
        final_results (dict): Same as input but with corrected umis.
        corrected_umis (int): How many umis have been corrected.
    """
    print('Correcting umis')
    corrected_umis = 0
    for cell_barcode in final_results:
        for TAG in final_results[cell_barcode]:
            if len(final_results[cell_barcode][TAG]) > 1:
                umi_clusters = network.UMIClusterer()
                UMIclusters = umi_clusters(
                    final_results[cell_barcode][TAG].keys(),
                    final_results[cell_barcode][TAG],
                    collapsing_threshold)
                for umi_cluster in UMIclusters:  # This is a list with the first element the dominant barcode
                    if(len(umi_cluster) > 1):  # This means we got a correction
                        major_umi = umi_cluster[0]
                        for minor_umi in umi_cluster[1:]:
                            corrected_umis += 1
                            temp = final_results[cell_barcode][TAG].pop(minor_umi)
                            final_results[cell_barcode][TAG][major_umi] += temp
    print('Corrected {} umis'.format(corrected_umis))
    return(final_results, corrected_umis)


def correct_cells(final_results, umis_per_cell, collapsing_threshold, expected_cells):
    """
    Corrects cell barcodes.
    
    Args:
        final_results (dict): Dict of dict of Counters with mapping results.
        umis_per_cell (Counter): Counter of number of umis per cell.
        collapsing_threshold (int): Max distance between umis.
        expected_cells (int): Number of expected cells.
    
    Returns:
        final_results (dict): Same as input but with corrected umis.
        corrected_umis (int): How many umis have been corrected.
    """
    print('Correcting cell barcodes')
    corrected_barcodes = 0
    try:
        cell_whitelist, true_to_false_map = umi_methods.getCellWhitelist(
            cell_barcode_counts=umis_per_cell,
            expect_cells=expected_cells,
            cell_number=False,
            error_correct_threshold=collapsing_threshold,
            plotfile_prefix=False)
        if true_to_false_map:
            for real_barcode in true_to_false_map:
                for fake_barcode in true_to_false_map[real_barcode]:
                    temp = final_results.pop(fake_barcode)
                    corrected_barcodes += 1
                    for TAG in temp.keys():
                        final_results[real_barcode][TAG].update(temp[TAG])
                    temp_umi_counts = umis_per_cell.pop(fake_barcode)
                    umis_per_cell[real_barcode] += temp_umi_counts
            print('Corrected {} cell barcodes'.format(corrected_barcodes))
    except Exception as e:
        print('Could not find a good local minima for correction.\nNo cell barcode correction was done.')
    return(final_results, umis_per_cell, corrected_barcodes)


def generate_sparse_matrices(final_results, ordered_tags_map, top_cells):
    """
    Create two sparse matrices with umi and read counts.

    Args:
        final_results (dict): Results in a dict of dicts of Counters.
        ordered_tags_map (dict): Tags in order with indexes as values.

    Returns:
        umi_results_matrix (scipy.sparse.dok_matrix): UMI counts
        read_results_matrix (scipy.sparse.dok_matrix): Read counts

    """
    umi_results_matrix = sparse.dok_matrix((len(ordered_tags_map) ,len(top_cells)), dtype=int16)
    read_results_matrix = sparse.dok_matrix((len(ordered_tags_map) ,len(top_cells)), dtype=int16)
    for i,cell_barcode in enumerate(top_cells):
        for j,TAG in enumerate(final_results[cell_barcode]):
            if final_results[cell_barcode][TAG]:
                umi_results_matrix[ordered_tags_map[TAG],i] = len(final_results[cell_barcode][TAG])
                read_results_matrix[ordered_tags_map[TAG],i] = sum(final_results[cell_barcode][TAG].values())
    return(umi_results_matrix, read_results_matrix)

