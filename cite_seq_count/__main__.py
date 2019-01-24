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
        prog='CITE-seq-Count', formatter_class=RawTextHelpFormatter,
        description=("This script counts matching antibody tags from two fastq "
                     "files. Version {}".format(version)),
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
    barcodes.add_argument('--umi_collapsing_dist', dest='umi_threshold',
                          required=False, type=int, default=2,
                          help="threshold for umi collapsing.")
    barcodes.add_argument('--bc_collapsing_dist', dest='bc_threshold',
                          required=False, type=int, default=1,
                          help="threshold for cellular barcode collapsing.")
    cells = parser.add_argument_group(
        'Cells',
        description=("Expected number of cells and potential whitelist")
    )

    cells.add_argument(
        '-cells', '--expected_cells', dest='expected_cells', required=True, type=int,
        help=("Number of expected cells from your run."), default=0
    )
    cells.add_argument(
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
        'TAG filters',
        description=("Filtering and trimming for read2.")
    )
    filters.add_argument(
        '--max-errors', dest='max_error',
        required=False, type=int, default=2,
        help=("Maximum Levenshtein distance allowed for antibody barcodes.")
    )
    
    filters.add_argument(
        '-trim', '--start-trim', dest='start_trim',
        required=False, type=int, default=0,
        help=("Number of bases to discard from read2.")
    )
    

    # Remaining arguments.
    parser.add_argument('-T', '--threads', required=False, type=int,
                        dest='n_threads', default=cpu_count(),
                        help="How many threads are to be used for running the program")
    parser.add_argument('-n', '--first_n', required=False, type=int,
                        dest='first_n', default=None,
                        help="Select N reads to run on instead of all.")
    parser.add_argument('-o', '--output', required=False, type=str, default='Results',
                        dest='outfolder', help="Results will be written to this folder")
    parser.add_argument('--dense', required=False, action='store_true', default=False,
                        dest='dense', help="Add a dense output to the results folder")
    parser.add_argument('-u', '--unmapped-tags', required=False, type=str,
                        dest='unmapped_file', default='unmapped.csv',
                        help="Write table of unknown TAGs to file.")
    parser.add_argument('-ut', '--unknown-top-tags', required=False,
                        dest='unknowns_top', type=int, default=100,
                        help="Top n unmapped TAGs.")
    parser.add_argument('--debug', action='store_true',
                        help="Print extra information for debugging.")
    parser.add_argument('--version', action='version', version='CITE-seq-Count v{}'.format(version),
                        help="Print version number.")
    
    # Finally! Too many options XD
    return parser


def create_report(n_reads, reads_per_cell, no_match, version, start_time, ordered_tags_map, umis_corrected, bcs_corrected, args):
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
    total_unmapped = sum(no_match.values())
    mapped_perc = round((total_mapped/n_reads)*100)
    unmapped_perc = round((total_unmapped/n_reads)*100)
    
    with open(os.path.join(args.outfolder, 'run_report.yaml'), 'w') as report_file:
        report_file.write(
"""Date: {}
Running time: {}
CITE-seq-Count Version: {}
Reads processed: {}
Percentage mapped: {}
Percentage unmapped: {}
Correction:
\tCell barcodes collapsing threshold: {}
\tCell barcodes corrected: {}
\tUMI collapsing threshold: {}
\tUMIs corrected: {}
Run parameters:
\tRead1_filename: {}
\tRead2_filename: {}
\tCell barcode:
\t\tFirst position: {}
\t\tLast position: {}
\tUMI barcode:
\t\tFirst position: {}
\t\tLast position: {}
\tExpected cells: {}
\tTags max errors: {}
\tStart trim: {}
""".format(
            datetime.datetime.today().strftime('%Y-%m-%d'),
            secondsToText.secondsToText(time.time()-start_time),
            version,
            n_reads,
            mapped_perc,
            unmapped_perc,
            args.bc_threshold,
            bcs_corrected,
            args.umi_threshold,
            umis_corrected,
            args.read1_path,
            args.read2_path,
            args.cb_first,
            args.cb_last,
            args.umi_first,
            args.umi_last,
            args.expected_cells,
            args.max_error,
            args.start_trim))

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
    ab_map = preprocessing.check_tags(ab_map, args.max_error)
    # Get reads length. So far, there is no validation for Read2.
    read1_length = preprocessing.get_read_length(args.read1_path)
    read2_length = preprocessing.get_read_length(args.read2_path)
    # Check Read1 length against CELL and UMI barcodes length.
    (
        barcode_slice, 
        umi_slice, 
        barcode_umi_length
    ) = preprocessing.check_barcodes_lengths(
            read1_length,
            args.cb_first,
            args.cb_last, 
            args.umi_first, args.umi_last)
    
    if args.first_n:
        n_lines = args.first_n*4
    else:
        n_lines = preprocessing.get_n_lines(args.read1_path)
    n_reads = int(n_lines/4)
    n_threads = args.n_threads
    
    print('Started mapping')
    #Run with one process
    if n_threads <= 1 or n_reads < 1000001:
        print('CITE-seq-Count is running with one core.')
        (
            final_results,
            merged_no_match) = processing.map_reads(
                read1_path=args.read1_path,
                read2_path=args.read2_path,               
                tags=ab_map,
                barcode_slice=barcode_slice,
                umi_slice=umi_slice,
                indexes=[0,n_reads],
                whitelist=whitelist,
                debug=args.debug,
                start_trim=args.start_trim,
                maximum_distance=args.max_error)
        print('Mapping done')
        umis_per_cell = Counter()
        reads_per_cell = Counter()
        for cell_barcode,counts in final_results.items():
            umis_per_cell[cell_barcode] = sum([len(counts[UMI]) for UMI in counts if UMI != 'unmapped'])
            reads_per_cell[cell_barcode] = sum([sum(counts[UMI].values()) for UMI in counts if UMI != 'unmapped'])
    else:
        # Run with multiple processes
        print('CITE-seq-Count is running with {} cores.'.format(n_threads))
        p = Pool(processes=n_threads)
        chunk_indexes = preprocessing.chunk_reads(n_reads, n_threads)
        parallel_results = []

        for indexes in chunk_indexes:
           p.apply_async(processing.map_reads,
                args=(
                    args.read1_path,
                    args.read2_path,
                    ab_map,
                    barcode_slice,
                    umi_slice,
                    indexes,
                    whitelist,
                    args.debug,
                    args.start_trim,
                    args.max_error),
                callback=parallel_results.append,
                error_callback=sys.stderr)
        p.close()
        p.join()
        print('Mapping done')
        print('Merging results')
        (
            final_results,
            umis_per_cell,
            reads_per_cell,
            merged_no_match
        ) = processing.merge_results(parallel_results=parallel_results)
        del(parallel_results)

    # Correct cell barcodes
    (
        final_results,
        umis_per_cell,
        bcs_corrected
    ) = processing.correct_cells(
            final_results=final_results,
            umis_per_cell=umis_per_cell,
            expected_cells=args.expected_cells,
            collapsing_threshold=args.bc_threshold)
    
    # Correct umi barcodes
    (
        final_results,
        umis_corrected
    ) = processing.correct_umis(
        final_results=final_results,
        collapsing_threshold=args.umi_threshold)

    ordered_tags_map = OrderedDict()
    for i,tag in enumerate(ab_map.values()):
        ordered_tags_map[tag] = i
    ordered_tags_map['unmapped'] = i + 1


    # Sort cells by number of mapped umis
    if not whitelist:
        top_cells_tuple = umis_per_cell.most_common(args.expected_cells)
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
    


    (
        umi_results_matrix,
        read_results_matrix
    ) = processing.generate_sparse_matrices(
        final_results=final_results,
        ordered_tags_map=ordered_tags_map,
        top_cells=top_cells)
    io.write_to_files(
        sparse_matrix=umi_results_matrix,
        top_cells=top_cells,
        ordered_tags_map=ordered_tags_map,
        data_type='umi',
        outfolder=args.outfolder)
    io.write_to_files(
        sparse_matrix=read_results_matrix,
        top_cells=top_cells,
        ordered_tags_map=ordered_tags_map,
        data_type='read',
        outfolder=args.outfolder)
      
    top_unmapped = merged_no_match.most_common(args.unknowns_top)
    with open(os.path.join(args.outfolder, args.unmapped_file),'w') as unknown_file:
        unknown_file.write('tag,count\n')
        for element in top_unmapped:
            unknown_file.write('{},{}\n'.format(element[0],element[1]))
    create_report(
        n_reads=n_reads,
        reads_per_cell=reads_per_cell,
        no_match=merged_no_match,
        version=version,
        start_time=start_time,
        ordered_tags_map=ordered_tags_map,
        umis_corrected=umis_corrected,
        bcs_corrected=bcs_corrected,
        args=args)
    if args.dense:
        print('Writing dense format output')
        io.write_dense(
            sparse_matrix=umi_results_matrix,
            index=list(ordered_tags_map.keys()),
            columns=top_cells,
            file_path=os.path.join(args.outfolder, 'dense_umis.tsv'))

if __name__ == '__main__':
    main()
