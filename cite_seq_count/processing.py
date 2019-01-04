import time
import gzip
import sys
import os
import Levenshtein
from collections import Counter
from itertools import islice

def classify_reads_multi_process(read1_path, read2_path, chunk_size,
                    tags, barcode_slice, umi_slice,
                    regex_pattern, first_line, whitelist,
                    legacy, debug):
    """Read through R1/R2 files and generate a islice starting at a specific index.

    It reads both Read1 and Read2 files, creating a dict based on cell barcode.

    Args:
        read_path1 (string): Path to R1.fastq.gz
        read_path2 (string): Path to R2.fastq.gz
        chunk_size (int): The number of lines to process 
        tags (dict): A dictionary with the TAGs + TAG Names.
        barcode_slice (slice): A slice for extracting the Barcode portion from the
            sequence.
        umi_slice (slice): A slice for extracting the UMI portion from the
            sequence.
        regex_pattern (regex.Pattern): An object that matches against any of the
            provided TAGs within the maximum distance provided.
        first_line (int): The first line to process in the file.
        whitelist (set): The set of white-listed barcodes.
        include_no_match (bool, optional): Whether to keep track of the
            `no_match` tags. Default is True.
        legacy (bool): If the tags contain a polyA or not.
        debug (bool): Print debug messages. Default is False.

    Returns:
        pandas.DataFrame: Matrix with the resulting counts.
        dict(int): A dictionary with the counts for each `no_match` TAG, based
            on the length of the longest provided TAG.

    """
    results_table = {}
    no_match_table = Counter()
    # Get the length of the longest TAG.
    longest_ab_tag = len(next(iter(tags)))
    n = 0
    t=time.time()
    if legacy:
        max_tag_length = longest_ab_tag + 6
    else:
        max_tag_length = longest_ab_tag
    with gzip.open(read1_path, 'rt') as textfile1, \
         gzip.open(read2_path, 'rt') as textfile2:
        
        # Read all 2nd lines from 4 line chunks. If first_n not None read only 4 times the given amount.
        secondlines = islice(zip(textfile1, textfile2), first_line, first_line + chunk_size - 1, 4)
        for read1, read2 in secondlines:
                n+=1
                if n % 1000000 == 0:
                    print("Processed 1,000,000 reads in {:.4} secondes. Total "
                        "reads: {:,} in child {}".format(time.time()-t, n, os.getpid()))
                    sys.stdout.flush()
                    t = time.time()
                read1 = read1.strip()
                read2 = read2.strip()
                cell_barcode = read1[barcode_slice]
                if whitelist is not None:
                    if cell_barcode not in whitelist:
                        continue
                if cell_barcode not in results_table:
                    results_table[cell_barcode] = {}
                    #for entry in entry_list:
                        #results_table[cell_barcode][entry]=Counter()
                TAG_seq = read2[:max_tag_length]
                UMI = read1[umi_slice]
                if debug:
                    print(
                        "\nline:{0}\n"
                        "cell_barcode:{1}\tUMI:{2}\tTAG_seq:{3}\n"
                        "line length:{4}\tcell barcode length:{5}\tUMI length:{6}\tTAG sequence length:{7}"
                        .format(read1 + read2, cell_barcode, UMI, TAG_seq,
                                len(read1 + read2), len(cell_barcode), len(UMI), len(TAG_seq)
                        )
                    )
                    sys.stdout.flush()
                    sys.stderr.flush()

                # Apply regex to Read2.
                match = regex_pattern.search(TAG_seq)
                if match:
                    # If a match is found, keep only the matching portion.
                    TAG_seq = match.group(0)
                    # Get the distance by adding up the errors found:
                    #   substitutions, insertions and deletions.
                    distance = sum(match.fuzzy_counts)
                    # To get the matching TAG, compare `match` against each TAG.
                    for tag, name in tags.items():
                        # This time, calculate the distance using the faster function
                        # `Levenshtein.distance` (which does the same). Thus, both
                        # determined distances should match.
                        if Levenshtein.distance(tag, TAG_seq) <= distance:
                            #results_table[cell_barcode]['total_reads'][UMI] += 1
                            try:
                                results_table[cell_barcode][name][UMI] += 1
                            except:
                                results_table[cell_barcode][name]=Counter()
                                results_table[cell_barcode][name][UMI] += 1
                            break
                
                else:
                    # No match
                    try:
                        results_table[cell_barcode]['no_match'][UMI] += 1
                    except:
                        results_table[cell_barcode]['no_match']=Counter()
                        results_table[cell_barcode]['no_match'][UMI] += 1
                    no_match_table[TAG_seq] += 1
    print("Counting done for process {}. Processed {:,} reads".format(os.getpid(), n))
    sys.stdout.flush()
    
    return(results_table, no_match_table)

def merge_results(parallel_results):
    """Merge chunked results from parallel processing.

    Args:
        parallel_results (list): List of dict with mapping results.

    Returns:
        merged_results (dict): Results combined as a dict of dict of Counters
        umis_per_cell (Counter): Total umis per cell as a Counter
        merged_no_match (Counter): Unmapped tags as a Counter
    """
    
    merged_results = {}
    
    merged_no_match = Counter()
    umis_per_cell = Counter()
    for chunk in parallel_results:
        mapped = chunk[0]
        unmapped = chunk[1]
        for cell_barcode in mapped:
            if cell_barcode not in merged_results:
                merged_results[cell_barcode] = dict()
            for TAG in mapped[cell_barcode]:
                # Test the counter. If empty, returns false
                if mapped[cell_barcode][TAG]:             
                    if TAG not in merged_results[cell_barcode]:
                        merged_results[cell_barcode][TAG] = Counter()
                    for UMI in mapped[cell_barcode][TAG]:
                        merged_results[cell_barcode][TAG][UMI] += mapped[cell_barcode][TAG][UMI]
                    umis_per_cell[cell_barcode] += len(mapped[cell_barcode][TAG])
        merged_no_match.update(unmapped)
    return(merged_results, umis_per_cell, merged_no_match)
