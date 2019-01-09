#!/usr/bin/env python3.6
"""
Author: Patrick Roelli
"""
import sys
import time
import os
import datetime
import pkg_resources

from argparse import ArgumentParser
from argparse import RawTextHelpFormatter
from collections import OrderedDict
from collections import Counter

from multiprocess import cpu_count
from multiprocess import Pool

from cite_seq_count import preprocessing
from cite_seq_count import processing
from cite_seq_count import io
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
    allowed_errors = ['i','s','d','e','i,s','i,d','s,d']
    filters.add_argument(
        '-e', '--error-type', dest='error_type',
        required=False, type=str, default='s',
        help=("Error type for the regex match."
            "\ni: Insertions only\ns: Substitutions only\n"
            "d: Deletion only\ne: All of the above"),
        choices=allowed_errors
    )
    
    # Remaining arguments.
    parser.add_argument('-T', '--threads', required=False, type=int,
                        dest='n_threads', default=cpu_count(),
                        help="How many threads are to be used for running the program")
    parser.add_argument('-n', '--first_n', required=False, type=int,
                        dest='first_n', default=None,
                        help="Select N reads to run on instead of all.")
    parser.add_argument('-o', '--output', required=True, type=str,
                        dest='outfolder', help="Results will be written to this folder")
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


def create_report(n_reads, reads_per_cell, no_match, version, start_time, ordered_tags_map, args):
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
    total_mapped = sum(reads_per_cell.values())
    print(total_mapped)
    total_unmapped = sum(no_match.values())
    mapped_perc = round((total_mapped/n_reads)*100)
    unmapped_perc = round((total_unmapped/n_reads)*100)
    
    with open(os.path.join(args.outfolder, 'run_report.yaml'), 'w') as report_file:
        report_file.write(
"""Date: {}
Running time: {}
CITE-seq-Count Version: {}
reads processed: {}
Percentage mapped: {}
Percentage unmapped: {}
Parameters:
\tRead1_filename: {}
\tRead2_filename: {}
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
        whitelist = preprocessing.parse_whitelist_csv(args.whitelist,
                                        args.cb_last - args.cb_first + 1)
    else:
        whitelist = None

    # Load TAGs/ABs.
    ab_map = preprocessing.parse_tags_csv(args.tags)
    ab_map = preprocessing.check_tags(ab_map, args.hamming_thresh)
    # Get reads length. So far, there is no validation for Read2.
    read1_length = preprocessing.get_read_length(args.read1_path)
    read2_length = preprocessing.get_read_length(args.read2_path)

    # Check Read1 length against CELL and UMI barcodes length.
    (barcode_slice, 
     umi_slice, 
     barcode_umi_length) = preprocessing.check_read_lengths(read1_length, args.cb_first,
                                              args.cb_last, 
                                              args.umi_first, args.umi_last)
    
    # Generate the compiled regex pattern.
    regex_pattern = preprocessing.generate_regex(ab_map, args.hamming_thresh, args.error_type, max_poly_a=6)

    
    if args.first_n:
        n_lines = args.first_n*4
    else:
        n_lines = preprocessing.get_n_lines(args.read1_path, args.first_n)  
    
    n_reads = n_lines/4
    n_threads = args.n_threads
    
    print('Started mapping')
    #Run with one process
    if n_threads <= 1:
        print('CITE-seq-Count is running with only one core.')
        (final_results, merged_no_match, total_reads) = processing.classify_reads_multi_process(
                args.read1_path,
                args.read2_path,               
                n_lines,
                ab_map,
                barcode_slice,
                umi_slice,
                regex_pattern,
                1,
                whitelist,
                args.legacy,
                args.debug)
        print('Mapping done')
        umis_per_cell = Counter()
        reads_per_cell = Counter()
        for cell_barcode,counts in final_results.items():
            umis_per_cell[cell_barcode] = sum([len(counts[UMI]) for UMI in counts if UMI != 'no_match' and UMI != 'bad_construct'])
            reads_per_cell[cell_barcode] = sum([sum(counts[UMI].values()) for UMI in counts if UMI != 'no_match' and UMI != 'bad_construct'])
    else:
        # Run with multiple processes
        p = Pool(processes=n_threads)
        # We need the chunk size to be a multiple of 4 to always start on the correct line in the fastq file
        chunk_size=round((n_lines-1)/n_threads)
        while chunk_size % 4 != 0:
            chunk_size +=1
        chunks = range(1,n_lines,chunk_size)

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
                callback=parallel_results.append,
                error_callback=sys.stderr)
        p.close()
        p.join()
        print('Mapping done')
        print('Merging results')
        (final_results, umis_per_cell, reads_per_cell, merged_no_match, total_reads) = processing.merge_results(parallel_results)                
        del(parallel_results)
    
    ordered_tags_map = OrderedDict()
    for i,tag in enumerate(ab_map.values()):
        ordered_tags_map[tag] = i
    ordered_tags_map['unmapped'] = i + 1


    # Sort cells by number of mapped umis
    if not whitelist:
        top_cells_tuple = umis_per_cell.most_common(args.cells)
        top_cells = set([pair[0] for pair in top_cells_tuple])
    else:
        top_cells = whitelist
        # Add potential missing cell barcodes.
        for missing_cell in whitelist:
            if missing_cell in final_results:
                continue
            else:
                final_results[missing_cell] = dict()
                for TAG in ordered_tags_map:
                    final_results[missing_cell][TAG] = 0
                top_cells.add(missing_cell)
    

    (umi_results_matrix, read_results_matrix) = processing.generate_sparse_matrices(final_results, ordered_tags_map, top_cells)
    
    io.write_to_files(umi_results_matrix, top_cells, ordered_tags_map, 'umi', args.outfolder)
    io.write_to_files(read_results_matrix, top_cells, ordered_tags_map, 'read', args.outfolder)
      
    # Save no_match TAGs to `args.unknowns_file` file.
    if args.unknowns_file:
        # Filter unknown TAGs base on the specified cutoff
        top_unmapped = merged_no_match.most_common(args.unknowns_top)
        with open(os.path.join(args.outfolder, args.unknowns_file),'w') as unknown_file:
            unknown_file.write('tag,count\n')
            for element in top_unmapped:
                unknown_file.write('{},{}\n'.format(element[0],element[1]))
    create_report(
        total_reads,
        reads_per_cell,
        merged_no_match,
        version,
        start_time,
        ordered_tags_map,
        args)

if __name__ == '__main__':
    main()
