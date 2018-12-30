#!/usr/bin/env python3
"""
Authors: Christoph Hafemeister, Patrick Roelli
"""
import csv
import gzip
import locale
import sys
import time
import warnings
import os

from argparse import ArgumentParser
from argparse import RawTextHelpFormatter
from collections import defaultdict
from collections import OrderedDict
from itertools import islice
from itertools import combinations

from multiprocess import Process, Manager, Pool
import Levenshtein
import pandas as pd
import pkg_resources
import regex
from scipy.sparse import dok_matrix
from scipy.io import mmwrite
from numpy import int16

version = pkg_resources.require("cite_seq_count")[0].version


def get_args():
    """
    Get args.
    """
    parser = ArgumentParser(
        prog='CITE Seq Count', formatter_class=RawTextHelpFormatter,
        description=("This script counts matching antibody tags from two fastq "
                     "files. Version {}".format(version))
    )

    # REQUIRED INPUTS group.
    inputs = parser.add_argument_group('Inputs',
                                       description="Required input files.")
    inputs.add_argument('-R1', '--read1', dest='read1_path', required=True,
                        help="The path of Read1 in gz format.")
    inputs.add_argument('-R2', '--read2', dest='read2_path', required=True,
                        help="The path of Read2 in gz format.")
    inputs.add_argument(
        '-t', '--tags', dest='tags', required=True,
        help=("The path to the csv file containing the antibody\n"
              "barcodes as well as their respective names.\n\n"
              "Example of an antibody barcode file structure:\n\n"
              "\tATGCGA,First_tag_name\n"
              "\tGTCATG,Second_tag_name")
    )

    # BARCODES group.
    barcodes = parser.add_argument_group(
        'Barcodes',
        description=("Positions of the cellular barcodes and UMI. If your "
                     "cellular barcodes and UMI\n are positioned as follows:\n"
                     "\tBarcodes from 1 to 16 and UMI from 17 to 26\n"
                     "then this is the input you need:\n"
                     "\t-cbf 1 -cbl 16 -umif 17 -umil 26")
    )
    barcodes.add_argument('-cbf', '--cell_barcode_first_base', dest='cb_first',
                          required=True, type=int,
                          help=("Postion of the first base of your cell "
                                "barcodes."))
    barcodes.add_argument('-cbl', '--cell_barcode_last_base', dest='cb_last',
                          required=True, type=int,
                          help=("Postion of the last base of your cell "
                                "barcodes."))
    barcodes.add_argument('-umif', '--umi_first_base', dest='umi_first',
                          required=True, type=int,
                          help="Postion of the first base of your UMI.")
    barcodes.add_argument('-umil', '--umi_last_base', dest='umi_last',
                          required=True, type=int,
                          help="Postion of the last base of your UMI.")

    # -cells and -whitelist are mutually exclusive options.
    barcodes_filtering = parser.add_mutually_exclusive_group(required=True)
    barcodes_filtering.add_argument(
        '-cells', '--expected_cells', dest='cells', required=False, type=int,
        help=("Number of expected cells from your run.")
    )
    barcodes_filtering.add_argument(
        '-wl', '--whitelist', dest='whitelist', required=False, type=str,
        help=("A csv file containning a whitelist of barcodes produced"
                      " by the mRNA data.\n\n"
                      "\tExample:\n"
                      "\tATGCTAGTGCTA\n\tGCTAGTCAGGAT\n\tCGACTGCTAACG\n\n"
                      "Or 10X-style:\n"
                      "\tATGCTAGTGCTA-1\n\tGCTAGTCAGGAT-1\n\tCGACTGCTAACG-1\n")
    )

    # FILTERS group.
    filters = parser.add_argument_group(
        'filters',
        description=("Filtering for structure of antibody barcodes as well as "
                    "maximum hamming\ndistance.")
    )
    filters.add_argument(
        '-hd', '--hamming-distance', dest='hamming_thresh',
        required=False, type=int, default=2,
        help=("Maximum hamming distance allowed for antibody barcode.")
    )
    
    # Remaining arguments.
    parser.add_argument('-n', '--first_n', required=False, type=int,
                        dest='first_n', default=None,
                        help="Select N reads to run on instead of all.")
    parser.add_argument('-o', '--output', required=True, type=str,
                        dest='outfile', help="Write result to file.")
    parser.add_argument('-u', '--unknown-tags', required=False, type=str,
                        dest='unknowns_file',
                        help="Write table of unknown TAGs to file.")
    parser.add_argument('-uc', '--unknown-tags-cutoff', required=False,
                        dest='unknowns_cutoff', type=int, default=10000,
                        help="Minimum counts to report an unknown TAG.")
    parser.add_argument('--debug', action='store_true',
                        help="Print extra information for debugging.")
    
    # REGEX related arguments.
    regex_pattern = parser.add_mutually_exclusive_group(required=False)
    regex_pattern.add_argument(
        '-tr', '--TAG_regex', dest='tag_regex', required=False, type=str,
        help=("Only use if you know what you are doing. The regex that will be "
              "used to validate an antibody barcode structure.\n"
              "Must be given in regex syntax.\n"
              "Example 1:\n"
              "\t\"^(GTCAACTCTTTAGCG|TGATGGCCTATTGGG)[TGC][A]{6,}\"\n"
              "\tMatches TAGs GTCAACTCTTTAGCG or TGATGGCCTATTGGG plus a T, G, "
              "or C, plus 6 or more As.\n"
              "Example 2:\n"
              "\"^[ATGC]{6}[TGC][A]{6,}\"\n"
              "Matches any 6 letter TAG.")
    )
    regex_pattern.add_argument(
        '-l', '--legacy', required=False, dest='legacy',
        default=False, action='store_true',
        help=("Use this option if you used an earlier version of the kit that "
              "adds a T,\nC, or G at the end of the sequence and you expect "
              "polyA tails in the data.")
    )

    # Finally! Too many options XD
    return parser


def parse_whitelist_csv(filename, barcode_length):
    """Reads white-listed barcodes from a CSV file.

    The function accepts plain barcodes or even 10X style barcodes with the
    `-1` at the end of each barcode.

    Args:
        filename (str): Whitelist barcode file.
        barcode_length (int): Length of the expected barcodes.

    Returns:
        set: The set of white-listed barcodes.

    """
    STRIP_CHARS = '0123456789- \t\n'
    with open(filename, mode='r') as csv_file:
        csv_reader = csv.reader(csv_file)
        whitelist = [row[0].strip(STRIP_CHARS) for row in csv_reader
                     if (len(row[0].strip(STRIP_CHARS)) == barcode_length)]
    return set(whitelist)


def parse_tags_csv(filename):
    """Reads the TAGs from a CSV file.

    The expected file format (no header) is: TAG,TAG_NAME.
    e.g. file content
        GTCAACTCTTTAGCG,Hashtag_1
        TGATGGCCTATTGGG,Hashtag_2
        TTCCGCCTCTCTTTG,Hashtag_3

    Args:
        filename (str): TAGs file.

    Returns:
        dict: A dictionary containing the TAGs and their names.

    """
    with open(filename, mode='r') as csv_file:
        csv_reader = csv.reader(csv_file)
        tags = {}
        for row in csv_reader:
            tags[row[0].strip()] = row[1].strip()
    return tags


def check_tags(tags, maximum_distance):
    """Evaluates the distance between the TAGs based on the `maximum distance`
    argument provided.

    Additionally, it adds the barcode to the name of the TAG circumventing the
    need of having to share the mapping of the antibody and the barcode.
    
    The output will have the keys sorted by TAG length (longer first). This
    way, longer barcodes will be evaluated first.

    Args:
        tags (dict): A dictionary with the TAGs + TAG Names.
        maximum_distance (int): The maximum Levenshtein distance allowed
            between two TAGs.

    Returns:
        collections.OrderedDict: An ordered dictionary containing the TAGs and
            their names in descendent order based on the length of the TAGs.

    """
    ordered_tags = OrderedDict()
    for tag in sorted(tags, key=len, reverse=True):
        ordered_tags[tag] = tags[tag] + '-' + tag

    # If only one TAG is provided, then no distances to compare.
    if (len(tags) == 1):
        return(ordered_tags)
    
    offending_pairs = []
    for a, b in combinations(ordered_tags.keys(), 2):
        distance = Levenshtein.distance(a, b)
        if (distance <= maximum_distance):
            offending_pairs.append([a, b, distance])
    
    # If offending pairs are found, print them all.
    if offending_pairs:
        print(
            '[ERROR] Minimum Levenshtein distance of TAGs barcode is less '
            'than given threshold.\n'
            'Please use a smaller distance.\n\n'
            'Offending case(s):\n'
        )
        for pair in offending_pairs:
            print(
                '\t{tag1}\n\t{tag2}\n\tDistance = {distance}\n'
                .format(tag1=pair[0], tag2=pair[1], distance=pair[2])
            )
        sys.exit('Exiting the application.\n')
    
    return(ordered_tags)


def get_read_length(filename):
    """Check wether SEQUENCE lengths are consistent in a FASTQ file and return
    the length.

    Args:
        filename (str): FASTQ file.

    Returns:
        int: The file's SEQUENCE length.

    """
    with gzip.open(filename, 'r') as fastq_file:
        secondlines = islice(fastq_file, 1, 1000, 4)
        temp_length = len(next(secondlines).rstrip())
        for sequence in secondlines:
            read_length = len(sequence.rstrip())
            if (temp_length != read_length):
                sys.exit(
                    '[ERROR] Sequence length is not consistent. Please, trim all '
                    'sequences at the same length.\n'
                    'Exiting the application.\n'
                )
    return(read_length)


def check_read_lengths(read1_length, cb_first, cb_last, umi_first, umi_last):
    """Check Read1 length against CELL and UMI barcodes length.

    Args:
        read1_length (int): Read1 length.
        cb_first (int): Barcode first base position for Read1.
        cb_last (int): Barcode last base position for Read1.
        umi_first (int): UMI first base position for Read1.
        umi_last (int): UMI last base position for Read1.

    Returns:
        slice: A `slice` object to extract the Barcode from the sequence string.
        slice: A `slice` object to extract the UMI from the sequence string.
        int: The Barcode + UMI length.

    """
    barcode_length = cb_last - cb_first + 1
    umi_length = umi_last - umi_first + 1
    barcode_umi_length = barcode_length + umi_length
    barcode_slice = slice(cb_first - 1, cb_last)
    umi_slice = slice(umi_first - 1, umi_last)

    if barcode_umi_length > read1_length:
        sys.exit(
            '[ERROR] Read1 length is shorter than the option you are using for '
            'Cell and UMI barcodes length. Please, check your options and rerun.\n\n'
            'Exiting the application.\n'
        )
    elif barcode_umi_length < read1_length:
        print(
            '[WARNING] Read1 length is {}bp but you are using {}bp for Cell '
            'and UMI barcodes combined.\nThis might lead to wrong cell '
            'attribution and skewed umi counts.\n'
            .format(read1_length, barcode_umi_length)
        )
    
    return(barcode_slice, umi_slice, barcode_umi_length)


def generate_regex(tags, maximum_distance, legacy=False, max_poly_a=6, read2_length=98, user_regex=None):
    """Generate regex based ont he provided TAGs.

    Args:
        tags (dict): A dictionary with the TAGs + TAG Names.
        maximum_distance (int): The maximum Levenshtein distance allowed
            between two TAGs.
        legacy (bool): `True` if you use an earlier version of the kit that adds
            a T, C, or G at the end and you expect polyA tails in the data.
            Default is False.
        max_poly_a (int): Run length of A's expected for the polyA tail. Default
            is 6.
        read2_length (int): Length of Read2. Default is 98.
        user_regex (str): A regular expression to use for TAG matching. Default
            is None.

    Returns:
        regex.Pattern: An object that matches against any of the provided TAGs
            within the maximum distance provided.

    """
    # Get a list of the available TAGs.
    tag_keys = tags.keys()

    # Get the length of the longest TAG.
    longest_ab_tag = len(next(iter(tags)))

    if user_regex:
        # If more than one TAG is provided and their length is different, issue a
        # warning.
        if len(tag_keys) > 1:
            for i in range(1, len(tag_keys)):
                if len(tag_keys[i]) != len(tag_keys[i - 1]):
                    print(
                        '[WARNING] Different length TAGs have been provided while '
                        'you specified a custom Regex. An OR method is recommended '
                        'for this scenarios. No additional validations will be '
                        'applied. Use it at your own risk.\n'
                    )
                    break
        
        regex_pattern = regex.compile(user_regex)
        return(regex_pattern)

    elif legacy:
        # Keep the minimum value between `max_poly_a` provided and the remaining
        # length of the read after removing the barcode length.
        polya_run = min(max_poly_a, read2_length - longest_ab_tag - 1)

        # Read comment below for `e` meaning in the regex.
        pattern = r'(^(\L<options>)[TGC][A]{{{},}}){{s<={}}}'.format(
            polya_run, maximum_distance
        )

    else:
        # `e` is part of the `regex` fuzzy logic: it means `error` in general,
        # whether it's a (s)ubstitution, (i)nsertion or (d)eletion. In this 
        # case, it means it allows `maximum_distance` errors to happen.
        pattern = r'(^\L<options>){{s<={}}}'.format(maximum_distance)

    # Compiling the regex makes it run faster.
    regex_pattern = regex.compile(pattern, options=tag_keys)
    return(regex_pattern)


def get_unique_lines(read1_filename, read2_filename, barcode_slice, umi_slice,
                     barcode_umi_length, first_n=None):
    """Read through R1/R2 files and generate a set without duplicate sequences.

    It reads both Read1 and Read2 files, creating a set based on Barcode + UMI
    + Read2 sequences. Note this means trimming Read1 after the UMI.

    Args:
        read1_filename (str): Read1 FASTQ file.
        read2_filename (str): Read2 FASTQ file.
        barcode_slice (slice): A slice for extracting the Barcode portion from the
            sequence.
        umi_slice (slice): A slice for extracting the UMI portion from the
            sequence.
        barcode_umi_length (int): The resulting length of adding the Barcode
            + UMI lengths.
        first_n (int, optional): The number of reads to subset from the FastQ
            files. Defaults to None.

    Returns:
        set: The unique combination of Barcode + UMI + Read2 sequences.

    """
    # Set object for storing unique Barcode+UMI+R2
    unique_lines = set()

    with gzip.open(read1_filename, 'rt') as textfile1, \
         gzip.open(read2_filename, 'rt') as textfile2:
        
        # Read all 2nd lines from 4 line chunks. If first_n not None read only 4 times the given amount.
        secondlines = islice(zip(textfile1, textfile2), 1, (first_n * 4 if first_n is not None else first_n), 4)
        print('loading')

        n = 0
        t = time.time()
        for read1, read2 in secondlines:
            read1 = read1.strip()
            read2 = read2.strip()
            line = read1[barcode_slice] + read1[umi_slice] + read2
            unique_lines.add(line)

            n += 1
            if n % 1000000 == 0:
                print("Loaded last 1,000,000 lines in {:.3} seconds. Total "
                      "lines loaded {:,} ".format(time.time()-t, n))
                t = time.time()

        print('{:,} reads loaded'.format(n))
        print('{:,} unique reads loaded'.format(len(unique_lines)))
    
    # Close R1/R2 files (release file handles)
    return(unique_lines)


def classify_reads(tags, unique_lines, barcode_slice, umi_slice,
                   barcode_umi_length, regex_pattern, whitelist=None,
                   include_no_match=True, debug=False):
    """Read through R1/R2 files and generate a set without duplicate sequences.

    It reads both Read1 and Read2 files, creating a set based on Barcode + UMI
    + Read2 sequences. Note this means trimming Read1 after the UMI.

    Args:
        tags (dict): A dictionary with the TAGs + TAG Names.
        unique_lines (set): The unique combination of Barcode + UMI + Read2
            sequences.
        barcode_slice (slice): A slice for extracting the Barcode portion from the
            sequence.
        umi_slice (slice): A slice for extracting the UMI portion from the
            sequence.
        barcode_umi_length (int): The resulting length of adding the Barcode
            + UMI lengths.
        regex_pattern (regex.Pattern): An object that matches against any of the
            provided TAGs within the maximum distance provided.
        whitelist (set): The set of white-listed barcodes.
        include_no_match (bool, optional): Whether to keep track of the
            `no_match` tags. Default is True.
        debug (bool): Print debug messages. Default is False.

    Returns:
        pandas.DataFrame: Matrix with the resulting counts.
        dict(int): A dictionary with the counts for each `no_match` TAG, based
            on the length of the longest provided TAG.

    """
    results_table = defaultdict(lambda: defaultdict(int))
    no_match_table = defaultdict(int)

    # Get the length of the longest TAG.
    longest_ab_tag = len(next(iter(tags)))

    n = 0
    t = time.time()
    for line in unique_lines:
        n += 1
        if n % 1000000 == 0:
            print("Processed 1,000,000 lines in {:.4} secondes. Total "
                  "lines processed: {:,}".format(time.time()-t, n))
            t = time.time()

        cell_barcode = line[barcode_slice]
        if whitelist:
            if cell_barcode not in whitelist:
                continue

        TAG_seq = line[barcode_umi_length:]
        if debug:
            UMI = line[umi_slice]
            print(
                "\nline:{0}\n"
                "cell_barcode:{1}\tUMI:{2}\tTAG_seq:{3}\n"
                "line length:{4}\tcell barcode length:{5}\tUMI length:{6}\tTAG sequence length:{7}"
                .format(line, cell_barcode, UMI, TAG_seq,
                        len(line), len(cell_barcode), len(UMI), len(TAG_seq)
                )
            )

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
                    results_table[cell_barcode]['total_reads'] += 1                    
                    results_table[cell_barcode][name][UMI] += 1
                    break
        
        else:
            # No match
            results_table[cell_barcode]['no_match'] += 1
            if include_no_match:
                tag = TAG_seq[:longest_ab_tag]
                no_match_table[tag] += 1
    
    print("Done counting")
    
    results_matrix = pd.DataFrame(results_table)
    if ('total_reads' not in results_matrix.index):
        exit('No match found. Please check your regex or tags file')
    
    return(results_matrix, no_match_table)

def info():
    print('{} reads were processed'.format(1000000))

def classify_reads_multi_process(chunk_id, tags, barcode_slice, umi_slice,
                   barcode_umi_length, regex_pattern, args, first_line, whitelist=None,
                   include_no_match=True, debug=False):
    """Read through R1/R2 files and generate a set without duplicate sequences.

    It reads both Read1 and Read2 files, creating a set based on Barcode + UMI
    + Read2 sequences. Note this means trimming Read1 after the UMI.

    Args:
        d (dict): A dictionnary for multiprocessing 
        tags (dict): A dictionary with the TAGs + TAG Names.
        unique_lines (set): The unique combination of Barcode + UMI + Read2
            sequences.
        barcode_slice (slice): A slice for extracting the Barcode portion from the
            sequence.
        umi_slice (slice): A slice for extracting the UMI portion from the
            sequence.
        barcode_umi_length (int): The resulting length of adding the Barcode
            + UMI lengths.
        regex_pattern (regex.Pattern): An object that matches against any of the
            provided TAGs within the maximum distance provided.
        whitelist (set): The set of white-listed barcodes.
        include_no_match (bool, optional): Whether to keep track of the
            `no_match` tags. Default is True.
        debug (bool): Print debug messages. Default is False.

    Returns:
        pandas.DataFrame: Matrix with the resulting counts.
        dict(int): A dictionary with the counts for each `no_match` TAG, based
            on the length of the longest provided TAG.

    """
    sys.stdout = open(str(os.getpid()) + ".out", "w")
    results_table = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
    no_match_table = defaultdict(int)
    # Get the length of the longest TAG.
    longest_ab_tag = len(next(iter(tags)))
    with gzip.open(args.read1_path, 'rt') as textfile1, \
         gzip.open(args.read2_path, 'rt') as textfile2:
        
        # Read all 2nd lines from 4 line chunks. If first_n not None read only 4 times the given amount.
        secondlines = islice(zip(textfile1, textfile2), first_line + 1, (args.first_n * 4 if args.first_n is not None else args.first_n), 4)
        for read1, read2 in secondlines:
                read1 = read1.strip()
                read2 = read2.strip()
                line = read1[barcode_slice] + read1[umi_slice] + read2


                cell_barcode = line[barcode_slice]
                if whitelist:
                    if cell_barcode not in whitelist:
                        continue

                TAG_seq = line[barcode_umi_length:]
                UMI = line[umi_slice]
                if debug:
                    print(
                        "\nline:{0}\n"
                        "cell_barcode:{1}\tUMI:{2}\tTAG_seq:{3}\n"
                        "line length:{4}\tcell barcode length:{5}\tUMI length:{6}\tTAG sequence length:{7}"
                        .format(line, cell_barcode, UMI, TAG_seq,
                                len(line), len(cell_barcode), len(UMI), len(TAG_seq)
                        )
                    )

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
                            results_table[cell_barcode]['total_reads'][UMI] += 1
                            results_table[cell_barcode][name][UMI] += 1
                            
                            break
                
                else:
                    # No match
                    results_table[cell_barcode]['no_match'][UMI] += 1
                    if include_no_match:
                        tag = TAG_seq[:longest_ab_tag]
                        no_match_table[tag] += 1
            
    # print("Done counting {} lines".format(len(unique_lines)))
    
    # results_matrix = pd.DataFrame(results_table)
    # if ('total_reads' not in results_matrix.index):
    #     exit('No match found. Please check your regex or tags file')
    #d[chunk_id]['results'] = results_table
    #d[chunk_id]['no_match'] = no_match_table
    
    # d['results{}'.format(chunk_id)] = results_table
    # d['no_match{}'.format(chunk_id)] = no_match_table
    #del(results_table)
    #del(no_match_table)
    info()
    #print(Processed reads from {} to {}'.format(i,i+1000000))
    return(results_table, no_match_table)

def merge_results(all_results_dict, ab_map):
    #final_results = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(int))))
    final_results = {}
    final_no_match = {}
    reads_per_cell = {}
    # chunk_ids = all_results_dict.keys()
    # chunk_res_map = dict()
    # for chunk_id in chunk_ids:
    #     chunk_res_map[chunk_id] = all_results_dict[chunk_id]['results'].keys()
    # all_results_dict[chunk_id]['results'].keys()
    for chunk in all_results_dict:
        mapped = chunk[0]
        unmapped = chunk[1]
        for cell_barcode in mapped:
            if not cell_barcode in final_results:
                final_results[cell_barcode] = dict()
                final_results[cell_barcode]['no_match'] = dict()
                final_results[cell_barcode]['total_reads'] = dict()
                reads_per_cell[cell_barcode] = 0
                for TAG in ab_map.values():
                    final_results[cell_barcode][TAG] = dict()            
            for TAG in mapped[cell_barcode]:
                for UMI in mapped[cell_barcode][TAG]:    
                    if UMI not in final_results[cell_barcode][TAG]:
                        final_results[cell_barcode][TAG][UMI]=0
                    final_results[cell_barcode][TAG][UMI] += mapped[cell_barcode][TAG][UMI]
                    reads_per_cell[cell_barcode] += mapped[cell_barcode][TAG][UMI]
    return(final_results, reads_per_cell)


def main():
    parser = get_args()
    if not sys.argv[1:]:
        parser.print_help(file=sys.stderr)
        sys.exit(2)

    # Parse arguments.
    args = parser.parse_args()
    if args.whitelist:
        whitelist = parse_whitelist_csv(args.whitelist,
                                        args.cb_last - args.cb_first + 1)
    else:
        whitelist = None

    # Load TAGs/ABs.
    ab_map = parse_tags_csv(args.tags)
    ab_map = check_tags(ab_map, args.hamming_thresh)
    
    # Get reads length. So far, there is no validation for Read2.
    read1_length = get_read_length(args.read1_path)
    #read2_length = get_read_length(args.read2_path)

    # Check Read1 length against CELL and UMI barcodes length.
    (barcode_slice, 
     umi_slice, 
     barcode_umi_length) = check_read_lengths(read1_length, args.cb_first,
                                              args.cb_last, 
                                              args.umi_first, args.umi_last)
    
    # Get unique combinations of Barcode+UMI+R2.
    # unique_lines = get_unique_lines(
    #     args.read1_path, args.read2_path, barcode_slice, umi_slice,
    #     barcode_umi_length, args.first_n)

    # Generate the compiled regex pattern.
    regex_pattern = generate_regex(ab_map, args.hamming_thresh, max_poly_a=6)

    #Initiate multiprocessing
    #manager = Manager()
    #d = manager.dict()
    nThreads = 4
    p = Pool(processes=nThreads)
    def blocks(files, size=65536):
        while True:
            b = files.read(size)
            if not b: break
            yield b
    
    if(args.first_n):
        nEntries = args.first_n
    else:
        print('Counting number of reads')
        with gzip.open(args.read1_path, "rt",encoding="utf-8",errors='ignore') as f:
            nEntries = int(sum(bl.count("\n") for bl in blocks(f))/4)
    chunks = range(0,nEntries,1000000)
    # Initiate proxy dicts for each chunk
    # for i in range(0,len(chunks)):
    #     d[i]=manager.dict()
    print('Started mapping')
    parallel_results = []
    for first_line,i in enumerate(range(0,len(chunks))):
        p.apply_async(classify_reads_multi_process,
            args=(i,ab_map, barcode_slice, umi_slice, barcode_umi_length, regex_pattern, args, i, whitelist, args.unknowns_file, args.debug),
            callback=parallel_results.append)
    p.close()
    p.join()
    print('Mapping done')
    # job = [
    # Process(
    #     target=classify_reads_multi_process,
    #     args=(
    #         d,
    #         i,
    #         ab_map,
    #         barcode_slice,
    #         umi_slice,
    #         barcode_umi_length,
    #         regex_pattern,
    #         args,
    #         first_line,
    #         whitelist,
    #         args.unknowns_file,
    #         args.debug)) for i,first_line in enumerate(chunks)]
    # _ = [p.start() for p in job]
    # _ = [p.join() for p in job]
    # # Perform the reads classification.
    # (results_matrix, no_match_table) = classify_reads_multi_process(
    #         ab_map, unique_lines, barcode_slice, umi_slice, barcode_umi_length,
    #         regex_pattern, whitelist, args.unknowns_file, args.debug)
    print('Merging results')
    (final_results,reads_per_cell) = merge_results(parallel_results, ab_map)
    del(parallel_results)

    if not whitelist:
        top_cells = sorted(reads_per_cell, key=reads_per_cell.get, reverse=True)[0:args.cells]
        final_results = {key: final_results[key] for key in final_results if key not in top_cells}
    else:
        #We just need to add the missing cell barcodes.
        for missing_cell in args.whitelist:
            if missing_cells in final_results:
                continue
            else:
                final_results[missing_cell] = dict()
                final_results[missing_cell]['no_match'] = dict()
                final_results[missing_cell]['total_reads'] = dict()
                reads_per_cell[missing_cell] = 0
                for TAG in ab_map.values():
                    final_results[missing_cell][TAG] = dict()            


    #Create an empty sparse matrix
    results_matrix = dok_matrix(((len(ab_map) + 2) ,len(final_results.keys())), dtype=int16)
    for i,cell_barcode in enumerate(final_results):
        for j,TAG in enumerate(final_results[cell_barcode]):
            value = len(final_results[cell_barcode][TAG])
            if value !=0:
                results_matrix[j,i]=value
    # # Add potential missing cells if whitelist is used.
    # if whitelist:
    #     results_matrix = results_matrix.reindex(whitelist, axis=1, fill_value=0)
    
    # # Replace any NA/NaN values with 0.
    # results_matrix.fillna(0, inplace=True)

    # # Keep only the TOP `args.cells` if provided.
    # if args.cells:
    #     most_reads_ordered = results_matrix.sort_values(by='total_reads',
    #                                                     ascending=False,
    #                                                     axis=1).columns
    #     n_top_cells = int(args.cells + args.cells/100*30)
    #     top_cells = most_reads_ordered[0:n_top_cells]
    #     results_matrix = results_matrix.loc[:, results_matrix.columns.isin(top_cells)]
    
    # Save results to `args.outfile` file.
    #results_matrix.to_csv(args.outfile, float_format='%.f')
    os.makedirs('mtx', exist_ok=True)
    mmwrite(os.path.join('mtx',args.outfile),results_matrix)
    with open(os.path.join('mtx','barcodes.csv'), 'w') as barcode_file:
        for barcode in final_results:
            barcode_file.write('{}\n'.format(barcode))
    with open(os.path.join('mtx','features.csv'), 'w') as feature_file:
        for barcode in final_results:
            for TAG in final_results[barcode]:
                feature_file.write('{}\n'.format(TAG))
            break

    
    # Save no_match TAGs to `args.unknowns_file` file.
    if args.unknowns_file:
        # Filter unknown TAGs base on the specified cutoff
        filtered_tags = {k:v
                         for k,v in no_match_table.items()
                         if v >= args.unknowns_cutoff}
        keys = list(filtered_tags.keys())
        vals = list(filtered_tags.values())
        no_match_matrix = pd.DataFrame({"tag": keys, "total": vals})
        no_match_matrix = no_match_matrix.sort_values(by='total', ascending=False)            
        no_match_matrix.to_csv(args.unknowns_file, float_format='%.f', index=False)


if __name__ == '__main__':
    main()
