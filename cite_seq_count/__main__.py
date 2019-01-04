#!/usr/bin/env python3
"""
Author: Patrick Roelli
"""
import csv
import gzip
import sys
import time
import os
import datetime

from argparse import ArgumentParser
from argparse import RawTextHelpFormatter
from collections import OrderedDict
from itertools import islice, combinations

import Levenshtein
from multiprocess import Process, cpu_count, Pool
import pkg_resources
import regex
from scipy import sparse
from scipy.io import mmwrite
from numpy import int16

from cite_seq_count import processing
from cite_seq_count import secondsToText

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
    parser.add_argument('-T', '--threads', required=False, type=int,
                        dest='n_threads', default=cpu_count(),
                        help="How many threads are to be used for running the program")
    parser.add_argument('-n', '--first_n', required=False, type=int,
                        dest='first_n', default=None,
                        help="Select N reads to run on instead of all.")
    parser.add_argument('-o', '--output', required=True, type=str,
                        dest='outfile', help="Write result to file.")
    parser.add_argument('-u', '--unknown-tags', required=False, type=str,
                        dest='unknowns_file',
                        help="Write table of unknown TAGs to file.")
    parser.add_argument('-ut', '--unknown-top-tags', required=False,
                        dest='unknowns_top', type=int, default=100,
                        help="Top n unmapped TAGs.")
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




def write_to_files(sparse_matrix, final_results, ordered_tags_map, data_type, outfile):
    """Write the sparse matrices to files in mtx format.

    Args:
        sparse_matrix (dok_matrix): Results in a sparse matrix.
        final_results (dict): Results in a dict of dicts of Counters.
        ordered_tags_map (dict): Tags in order with indexes as values.
        data_type (string): A string definning if the data is umi or read based.
        oufile (string): Path to the mtx file.
    
    """
    prefix = data_type + '_count'
    os.makedirs(prefix, exist_ok=True)
    mmwrite(os.path.join(prefix,outfile),sparse_matrix)
    with open(os.path.join(prefix,'barcodes.csv'), 'w') as barcode_file:
        for barcode in final_results:
            barcode_file.write('{}\n'.format(barcode))
    with open(os.path.join(prefix,'features.csv'), 'w') as feature_file:
        for feature in ordered_tags_map:
            feature_file.write('{}\n'.format(feature))
    

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
                umi_results_matrix[ordered_tags_map[TAG],i]=len(final_results[cell_barcode][TAG])
                read_results_matrix[ordered_tags_map[TAG],i]=sum(final_results[cell_barcode][TAG].values())
    return(umi_results_matrix, read_results_matrix)


def blocks(files, size=65536):
    """
    A fast way of counting the lines of a large file.
    Ref:
        https://stackoverflow.com/a/9631635/9178565

    Args:
        files (io.handler): A file handler 
        size (int): Block size
    Returns:
        A generator
    """
    while True:
        b = files.read(size)
        if not b: break
        yield b


def get_n_lines(file_path, top_n):
    """
    Determines how many lines have to be processed
    depending on options and number of available lines.
    Checks that the number of lines is a multiple of 4.

    Args:
        file_path (string): Path to a fastq.gz file
        top_n (int): Number of reads to be used as defined by the user.

    Returns:
        n_lines (int): Number of lines to be used
    """
    print('Counting number of reads')
    with gzip.open(file_path, "rt",encoding="utf-8",errors='ignore') as f:
        n_lines = sum(bl.count("\n") for bl in blocks(f))
    if n_lines %4 !=0:
        sys.exit('{}\'s number of lines is not a multiple of 4. The file might be corrupted.\n Exiting')
    if top_n and top_n*4 < n_lines:
        n_lines = top_n*4
    elif top_n and top_n*4 > n_lines:
        print('The data only contains {} reads. Processing all the reads'.format(n_lines/4))
    return(n_lines)


def create_report(n_reads, reads_matrix, no_match, version, start_time, ordered_tags_map, args):
    """
    Creates a report with details about the run in a yaml format.

    Args:
        n_reads (int): Number of reads that have been processed.
        reads_matrix (scipy.sparse.dok_matrix): A sparse matrix continining read counts.
        no_match (Counter): Counter of unmapped tags.
        version (string): CITE-seq-Count package version.
        start_time (time): Start time of the run.
        args (arg_parse): Arguments provided by the user.

    """
    ordered_tags_map.pop('no_match')
    total_mapped = 0
    for i in range(1,len(ordered_tags_map)):
        total_mapped += sum(reads_matrix[i,].values())
    total_unmapped = sum(no_match.values())
    mapped_perc = round((total_mapped/n_reads)*100)
    unmapped_perc = round((total_unmapped/n_reads)*100)
    
    with open('run_report.yaml', 'w') as report_file:
        report_file.write(
"""Date: {}
Running time: {}
CITE-seq-Count Version: {}
reads processed: {}
Percentage mapped: {}
Percentage unmapped: {}
Parameters:
\tRead1: {}
\tRead2: {}
\tCell barcode:
\t\tFirst position: {}
\t\tLast position: {}
\tUMI barcode:
\t\tFirst position: {}
\t\tLast position: {}
\tHamming distance: {}
\tLegacy: {}        
""".format(
            datetime.datetime.today().strftime('%Y-%m-%d'),
            secondsToText.secondsToText(time.time()-start_time),
            version,
            n_reads,
            mapped_perc,
            unmapped_perc,
            args.read1_path,
            args.read2_path,
            args.cb_first,
            args.cb_last,
            args.umi_first,
            args.umi_last,
            args.hamming_thresh,
            args.legacy))

def main():
    start_time = time.time()
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

    
    n_lines = get_n_lines(args.read1_path, args.first_n)  
    
    n_reads = n_lines/4
    n_threads = args.n_threads
    #Initiate multiprocessing leaving one thread available.
    p = Pool(processes=n_threads)
    # We need the chunk size to be a multiple of 4 to always start of the right line in the fastq file
    chunk_size=round((n_lines-1)/n_threads)
    while chunk_size % 4 != 0:
        chunk_size +=1
    chunks = range(1,n_lines,chunk_size)
    print('Started mapping')
    parallel_results = []
    for first_line in list(chunks):
       p.apply_async(processing.classify_reads_multi_process,
            args=(
                args.read1_path,
                args.read2_path,               
                chunk_size,
                ab_map,
                barcode_slice,
                umi_slice,
                regex_pattern,
                first_line,
                whitelist,
                args.legacy,
                args.debug),
            callback=parallel_results.append)
    p.close()
    p.join()

    print('Mapping done')

    print('Merging results')
    (final_results,reads_per_cell, merged_no_match) = processing.merge_results(parallel_results)
    del(parallel_results)
    
    ordered_tags_map = OrderedDict()
    for i,tag in enumerate(ab_map.values()):
        ordered_tags_map[tag] = i
    ordered_tags_map['no_match'] = i + 1

    # if not whitelist:
    top_cells_tuple = reads_per_cell.most_common(args.cells)
    top_cells = [pair[0] for pair in top_cells_tuple]
    #     final_results = {key: final_results[key] for key in top_cells}

    #We just need to add the missing cell barcodes.
    if whitelist:
        for missing_cell in args.whitelist:
            if missing_cell in final_results:
                continue
            else:
                final_results[missing_cell] = dict()
                for TAG in ordered_tags_map:
                    final_results[missing_cell][TAG] = 0
                top_cells.append(missing_cell)

    (umi_results_matrix, read_results_matrix) = generate_sparse_matrices(final_results, ordered_tags_map, top_cells)
    
    write_to_files(umi_results_matrix, final_results, ordered_tags_map, 'umi', args.outfile)
    write_to_files(read_results_matrix, final_results, ordered_tags_map, 'read', args.outfile)
      
    # Save no_match TAGs to `args.unknowns_file` file.
    if args.unknowns_file:
        # Filter unknown TAGs base on the specified cutoff
        top_unmapped = merged_no_match.most_common(args.unknowns_top)
        with open(args.unknowns_file,'w') as unknown_file:
            unknown_file.write('tag,count\n')
            for element in top_unmapped:
                unknown_file.write('{},{}\n'.format(element[0],element[1]))
    create_report(
        int(n_lines/4),
        read_results_matrix,
        merged_no_match,
        version,
        start_time,
        ordered_tags_map,
        args)

if __name__ == '__main__':
    main()
