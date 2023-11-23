#!/usr/bin/env python3.6
"""
Author: Patrick Roelli
"""
import sys
import time
import os
import datetime
import pkg_resources
import logging

from argparse import ArgumentParser
from argparse import RawTextHelpFormatter
from collections import OrderedDict
from collections import Counter
from collections import defaultdict

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
        prog="CITE-seq-Count",
        formatter_class=RawTextHelpFormatter,
        description=(
            "This script counts matching antibody tags from paired fastq "
            "files. Version {}".format(version)
        ),
    )

    # REQUIRED INPUTS group.
    inputs = parser.add_argument_group("Inputs", description="Required input files.")
    inputs.add_argument(
        "-R1",
        "--read1",
        dest="read1_path",
        required=True,
        help=(
            "The path of Read1 in gz format, or a comma-separated list of paths to all Read1 files in"
            " gz format (E.g. A1.fq.gz,B1.fq,gz,..."
        ),
    )
    inputs.add_argument(
        "-R2",
        "--read2",
        dest="read2_path",
        required=True,
        help=(
            "The path of Read2 in gz format, or a comma-separated list of paths to all Read2 files in"
            " gz format (E.g. A2.fq.gz,B2.fq,gz,..."
        ),
    )
    inputs.add_argument(
        "-t",
        "--tags",
        dest="tags",
        required=True,
        help=(
            "The path to the csv file containing the antibody\n"
            "barcodes as well as their respective names.\n\n"
            "Example of an antibody barcode file structure:\n\n"
            "\tATGCGA,First_tag_name\n"
            "\tGTCATG,Second_tag_name"
        ),
    )

    # BARCODES group.
    barcodes = parser.add_argument_group(
        "Barcodes",
        description=(
            "Positions of the cellular barcodes and UMI. If your "
            "cellular barcodes and UMI\n are positioned as follows:\n"
            "\tBarcodes from 1 to 16 and UMI from 17 to 26\n"
            "then this is the input you need:\n"
            "\t-cbf 1 -cbl 16 -umif 17 -umil 26"
        ),
    )
    barcodes.add_argument(
        "-cbf",
        "--cell_barcode_first_base",
        dest="cb_first",
        required=True,
        type=int,
        help=("Postion of the first base of your cell " "barcodes."),
    )
    barcodes.add_argument(
        "-cbl",
        "--cell_barcode_last_base",
        dest="cb_last",
        required=True,
        type=int,
        help=("Postion of the last base of your cell " "barcodes."),
    )
    barcodes.add_argument(
        "-umif",
        "--umi_first_base",
        dest="umi_first",
        required=True,
        type=int,
        help="Postion of the first base of your UMI.",
    )
    barcodes.add_argument(
        "-umil",
        "--umi_last_base",
        dest="umi_last",
        required=True,
        type=int,
        help="Postion of the last base of your UMI.",
    )
    barcodes.add_argument(
        "--umi_collapsing_dist",
        dest="umi_threshold",
        required=False,
        type=int,
        default=2,
        help="threshold for umi collapsing.",
    )
    barcodes.add_argument(
        "--no_umi_correction",
        required=False,
        action="store_true",
        default=False,
        dest="no_umi_correction",
        help="Deactivate UMI collapsing",
    )
    barcodes.add_argument(
        "--bc_collapsing_dist",
        dest="bc_threshold",
        required=False,
        type=int,
        default=1,
        help="threshold for cellular barcode collapsing.",
    )
    cells = parser.add_argument_group(
        "Cells", description=("Expected number of cells and potential whitelist")
    )

    cells.add_argument(
        "-cells",
        "--expected_cells",
        dest="expected_cells",
        required=True,
        type=int,
        help=("Number of expected cells from your run."),
        default=0,
    )
    cells.add_argument(
        "-wl",
        "--whitelist",
        dest="whitelist",
        required=False,
        type=str,
        help=(
            "A csv file containning a whitelist of barcodes produced"
            " by the mRNA data.\n\n"
            "\tExample:\n"
            "\tATGCTAGTGCTA\n\tGCTAGTCAGGAT\n\tCGACTGCTAACG\n\n"
            "Or 10X-style:\n"
            "\tATGCTAGTGCTA-1\n\tGCTAGTCAGGAT-1\n\tCGACTGCTAACG-1\n"
        ),
    )

    # FILTERS group.
    filters = parser.add_argument_group(
        "TAG filters", description=("Filtering and trimming for read2.")
    )
    filters.add_argument(
        "--max-errors",
        dest="max_error",
        required=False,
        type=int,
        default=2,
        help=("Maximum Levenshtein distance allowed for antibody barcodes."),
    )

    filters.add_argument(
        "-trim",
        "--start-trim",
        dest="start_trim",
        required=False,
        type=int,
        default=0,
        help=("Number of bases to discard from read2."),
    )

    filters.add_argument(
        "--sliding-window",
        dest="sliding_window",
        required=False,
        default=False,
        action="store_true",
        help=("Allow for a sliding window when aligning."),
    )

    # Read ID output
    readids = parser.add_argument_group(
        "Read IDs", description=("parameters to get read ids per cell")
    )

    readids.add_argument(
        "--store-read-ids",
        dest="store_read_ids",
        required=False,
        default=False,
        action="store_true",
        help=("Write read ids per cell to file."),
    )

    readids.add_argument(
        "--read-ids-whitelist",
        dest="read_ids_whitelist",
        required=False,
        type=str,
        help=(
            "A csv file containning a whitelist of barcodes for which"
            " read ids should be written to files.\n\n"
            "\tExample:\n"
            "\tATGCTAGTGCTA\n\tGCTAGTCAGGAT\n\tCGACTGCTAACG\n\n"
            "Or 10X-style:\n"
            "\tATGCTAGTGCTA-1\n\tGCTAGTCAGGAT-1\n\tCGACTGCTAACG-1\n"
        ),
    )

    # Remaining arguments.
    parser.add_argument(
        "-T",
        "--threads",
        required=False,
        type=int,
        dest="n_threads",
        default=cpu_count(),
        help="How many threads are to be used for running the program",
    )
    parser.add_argument(
        "-n",
        "--first_n",
        required=False,
        type=int,
        dest="first_n",
        default=None,
        help="Select N reads to run on instead of all.",
    )
    parser.add_argument(
        "-o",
        "--output",
        required=False,
        type=str,
        default="Results",
        dest="outfolder",
        help="Results will be written to this folder",
    )
    parser.add_argument(
        "--dense",
        required=False,
        action="store_true",
        default=False,
        dest="dense",
        help="Add a dense output to the results folder",
    )
    parser.add_argument(
        "-u",
        "--unmapped-tags",
        required=False,
        type=str,
        dest="unmapped_file",
        default="unmapped.csv",
        help="Write table of unknown TAGs to file.",
    )
    parser.add_argument(
        "-ut",
        "--unknown-top-tags",
        required=False,
        dest="unknowns_top",
        type=int,
        default=100,
        help="Top n unmapped TAGs.",
    )
    parser.add_argument(
        "--debug", action="store_true", help="Print extra information for debugging."
    )
    parser.add_argument(
        "--version",
        action="version",
        version="CITE-seq-Count v{}".format(version),
        help="Print version number.",
    )
    # Finally! Too many options XD
    return parser


def create_report(
    n_reads,
    reads_per_cell,
    no_match,
    version,
    start_time,
    ordered_tags_map,
    umis_corrected,
    bcs_corrected,
    bad_cells,
    args,
):
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
    total_unmapped = sum(no_match.values())
    total_mapped = sum(reads_per_cell.values()) - total_unmapped
    mapped_perc = round((total_mapped / n_reads) * 100)
    unmapped_perc = round((total_unmapped / n_reads) * 100)

    with open(os.path.join(args.outfolder, "run_report.yaml"), "w") as report_file:
        report_file.write(
            """Date: {}
Running time: {}
CITE-seq-Count Version: {}
Reads processed: {}
Percentage mapped: {}
Percentage unmapped: {}
Uncorrected cells: {}
Correction:
\tCell barcodes collapsing threshold: {}
\tCell barcodes corrected: {}
\tUMI collapsing threshold: {}
\tUMIs corrected: {}
Run parameters:
\tRead1_paths: {}
\tRead2_paths: {}
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
                datetime.datetime.today().strftime("%Y-%m-%d"),
                secondsToText.secondsToText(time.time() - start_time),
                version,
                n_reads,
                mapped_perc,
                unmapped_perc,
                len(bad_cells),
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
                args.start_trim,
            )
        )


def main():
    # Create logger and stream handler
    logger = logging.getLogger("cite_seq_count")
    logger.setLevel(logging.CRITICAL)
    ch = logging.StreamHandler()
    ch.setLevel(logging.CRITICAL)
    formatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )
    ch.setFormatter(formatter)
    logger.addHandler(ch)

    start_time = time.time()
    parser = get_args()
    if not sys.argv[1:]:
        parser.print_help(file=sys.stderr)
        sys.exit(2)

    # Parse arguments.
    args = parser.parse_args()
    if args.whitelist:
        print("Loading whitelist")
        (whitelist, args.bc_threshold) = preprocessing.parse_whitelist_csv(
            filename=args.whitelist,
            barcode_length=args.cb_last - args.cb_first + 1,
            collapsing_threshold=args.bc_threshold,
        )
    else:
        whitelist = False

    if args.read_ids_whitelist and args.store_read_ids:
        (read_ids_whitelist, _) = preprocessing.parse_whitelist_csv(
            filename=args.read_ids_whitelist,
            barcode_length=args.cb_last - args.cb_first + 1,
            collapsing_threshold=args.bc_threshold,
        )
    else:
        read_ids_whitelist = False

    # Load TAGs/ABs.
    ab_map = preprocessing.parse_tags_csv(args.tags)
    ab_map = preprocessing.check_tags(ab_map, args.max_error)

    # Identify input file(s)
    read1_paths, read2_paths = preprocessing.get_read_paths(
        args.read1_path, args.read2_path
    )

    # preprocessing and processing occur in separate loops so the program can crash earlier if
    # one of the inputs is not valid.
    read1_lengths = []
    read2_lengths = []
    for read1_path, read2_path in zip(read1_paths, read2_paths):
        # Get reads length. So far, there is no validation for Read2.
        read1_lengths.append(preprocessing.get_read_length(read1_path))
        read2_lengths.append(preprocessing.get_read_length(read2_path))
        # Check Read1 length against CELL and UMI barcodes length.
        (
            barcode_slice,
            umi_slice,
            barcode_umi_length,
        ) = preprocessing.check_barcodes_lengths(
            read1_lengths[-1],
            args.cb_first,
            args.cb_last,
            args.umi_first,
            args.umi_last,
        )
    # Ensure all files have the same input length
    # if len(set(read1_lengths)) != 1:
    #    sys.exit('Input barcode fastqs (read1) do not all have same length.\nExiting')

    # Initialize the counts dicts that will be generated from each input fastq pair
    final_results = defaultdict(lambda: defaultdict(Counter))
    readsid_per_cell = {}
    umis_per_cell = Counter()
    reads_per_cell = Counter()
    merged_no_match = Counter()
    number_of_samples = len(read1_paths)
    n_reads = 0

    # Print a statement if multiple files are run.
    if number_of_samples != 1:
        print("Detected {} files to run on.".format(number_of_samples))

    for read1_path, read2_path in zip(read1_paths, read2_paths):
        if args.first_n:
            n_lines = (args.first_n * 4) / number_of_samples
        else:
            n_lines = preprocessing.get_n_lines(read1_path)
        n_reads += int(n_lines / 4)
        n_threads = args.n_threads
        print("Started mapping")
        print("Processing {:,} reads".format(n_reads))
        # Run with one process
        if n_threads <= 1 or n_reads < 1000001:
            print("CITE-seq-Count is running with one core.")
            (_final_results, _readsid_per_cell, _merged_no_match) = processing.map_reads(
                read1_path=read1_path,
                read2_path=read2_path,
                tags=ab_map,
                barcode_slice=barcode_slice,
                umi_slice=umi_slice,
                indexes=[0, n_reads],
                whitelist=whitelist,
                debug=args.debug,
                start_trim=args.start_trim,
                maximum_distance=args.max_error,
                sliding_window=args.sliding_window,
                store_read_ids = args.store_read_ids,
            )
            print("Mapping done")
            _umis_per_cell = Counter()
            _reads_per_cell = Counter()
            for cell_barcode, counts in _final_results.items():
                _umis_per_cell[cell_barcode] = sum([len(counts[UMI]) for UMI in counts])
                _reads_per_cell[cell_barcode] = sum(
                    [sum(counts[UMI].values()) for UMI in counts]
                )
        else:
            # Run with multiple processes
            print("CITE-seq-Count is running with {} cores.".format(n_threads))
            p = Pool(processes=n_threads)
            chunk_indexes = preprocessing.chunk_reads(n_reads, n_threads)
            parallel_results = []

            for indexes in chunk_indexes:
                p.apply_async(
                    processing.map_reads,
                    args=(
                        read1_path,
                        read2_path,
                        ab_map,
                        barcode_slice,
                        umi_slice,
                        indexes,
                        whitelist,
                        args.debug,
                        args.start_trim,
                        args.max_error,
                        args.sliding_window,
                        args.store_read_ids,
                    ),
                    callback=parallel_results.append,
                    error_callback=sys.stderr,
                )
            p.close()
            p.join()
            print("Mapping done")
            print("Merging results")

            (
                _final_results,
                _readsid_per_cell,
                _umis_per_cell,
                _reads_per_cell,
                _merged_no_match,
            ) = processing.merge_results(parallel_results=parallel_results)
            del parallel_results

        # Update the overall counts dicts
        umis_per_cell.update(_umis_per_cell)
        reads_per_cell.update(_reads_per_cell)
        merged_no_match.update(_merged_no_match)
        for cell_barcode in _final_results:
            for tag in _final_results[cell_barcode]:
                if tag in final_results[cell_barcode]:
                    # Counter + Counter = Counter
                    final_results[cell_barcode][tag] += _final_results[cell_barcode][
                        tag
                    ]
                else:
                    # Explicitly save the counter to that tag
                    final_results[cell_barcode][tag] = _final_results[cell_barcode][tag]
        for cell_barcode in _readsid_per_cell:
            if cell_barcode not in readsid_per_cell:
                readsid_per_cell[cell_barcode] = []
            readsid_per_cell[cell_barcode] += _readsid_per_cell[cell_barcode]
    ordered_tags_map = OrderedDict()
    for i, tag in enumerate(ab_map.values()):
        ordered_tags_map[tag] = i
    ordered_tags_map["unmapped"] = i + 1

    # Correct cell barcodes
    if args.bc_threshold > 0:
        if len(umis_per_cell) <= args.expected_cells:
            print(
                "Number of expected cells, {}, is higher "
                "than number of cells found {}.\nNot performing"
                "cell barcode correction"
                "".format(args.expected_cells, len(umis_per_cell))
            )
            bcs_corrected = 0
        else:
            print("Correcting cell barcodes")
            if not whitelist:
                (
                    final_results,
                    final_readsid,
                    umis_per_cell,
                    bcs_corrected,
                ) = processing.correct_cells(
                    final_results=final_results,
                    readsid_per_cell=readsid_per_cell,
                    reads_per_cell=reads_per_cell,
                    umis_per_cell=umis_per_cell,
                    expected_cells=args.expected_cells,
                    collapsing_threshold=args.bc_threshold,
                    ab_map=ordered_tags_map,
                )
            else:
                (
                    final_results,
                    final_readsid,
                    umis_per_cell,
                    bcs_corrected,
                ) = processing.correct_cells_whitelist(
                    final_results=final_results,
                    readsid_per_cell=readsid_per_cell,
                    umis_per_cell=umis_per_cell,
                    whitelist=whitelist,
                    collapsing_threshold=args.bc_threshold,
                    ab_map=ordered_tags_map,
                )
    else:
        bcs_corrected = 0

    # If given, use whitelist for top cells
    if whitelist:
        top_cells = whitelist
        # Add potential missing cell barcodes.
        for missing_cell in whitelist:
            if missing_cell in final_results:
                continue
            else:
                final_results[missing_cell] = dict()
                for TAG in ordered_tags_map:
                    final_results[missing_cell][TAG] = Counter()
                top_cells.add(missing_cell)
    else:
        # Select top cells based on total umis per cell
        top_cells_tuple = umis_per_cell.most_common(args.expected_cells)
        top_cells = set([pair[0] for pair in top_cells_tuple])

    # UMI correction

    if args.no_umi_correction:
        # Don't correct
        umis_corrected = 0
        aberrant_cells = set()
    else:
        # Correct UMIS
        (final_results, umis_corrected, aberrant_cells) = processing.correct_umis(
            final_results=final_results,
            collapsing_threshold=args.umi_threshold,
            top_cells=top_cells,
            max_umis=20000,
        )

    # Remove aberrant cells from the top cells
    for cell_barcode in aberrant_cells:
        top_cells.remove(cell_barcode)

    # Create sparse aberrant cells matrix
    (umi_aberrant_matrix, read_aberrant_matrix) = processing.generate_sparse_matrices(
        final_results=final_results,
        ordered_tags_map=ordered_tags_map,
        top_cells=aberrant_cells,
    )

    # Write uncorrected cells to dense output
    if (len(aberrant_cells) > 0):
        io.write_dense(
            sparse_matrix=umi_aberrant_matrix,
            index=list(ordered_tags_map.keys()),
            columns=aberrant_cells,
            outfolder=os.path.join(args.outfolder, "uncorrected_cells"),
            filename="dense_umis.tsv",
        )

    # Create sparse matrices for results
    (umi_results_matrix, read_results_matrix) = processing.generate_sparse_matrices(
        final_results=final_results,
        ordered_tags_map=ordered_tags_map,
        top_cells=top_cells,
    )

    # Write umis to file
    io.write_to_files(
        sparse_matrix=umi_results_matrix,
        top_cells=top_cells,
        ordered_tags_map=ordered_tags_map,
        data_type="umi",
        outfolder=args.outfolder,
    )

    # Write reads to file
    io.write_to_files(
        sparse_matrix=read_results_matrix,
        top_cells=top_cells,
        ordered_tags_map=ordered_tags_map,
        data_type="read",
        outfolder=args.outfolder,
    )

    # Write unmapped sequences
    io.write_unmapped(
        merged_no_match=merged_no_match,
        top_unknowns=args.unknowns_top,
        outfolder=args.outfolder,
        filename=args.unmapped_file,
    )

    # Create report and write it to disk
    create_report(
        n_reads=n_reads,
        reads_per_cell=reads_per_cell,
        no_match=merged_no_match,
        version=version,
        start_time=start_time,
        ordered_tags_map=ordered_tags_map,
        umis_corrected=umis_corrected,
        bcs_corrected=bcs_corrected,
        bad_cells=aberrant_cells,
        args=args,
    )

    # Write dense matrix to disk if requested
    if args.dense:
        print("Writing dense format output")
        io.write_dense(
            sparse_matrix=umi_results_matrix,
            index=list(ordered_tags_map.keys()),
            columns=top_cells,
            outfolder=args.outfolder,
            filename="dense_umis.tsv",
        )
    
    # Write read ids to disk if requested
    if args.store_read_ids:
        print("Writting read ids")
        io.write_read_ids(
            readsid_per_cell = final_readsid,
            read_ids_whitelist = read_ids_whitelist,
            outfolder=args.outfolder,
            subfolder="read_ids"
        )


if __name__ == "__main__":
    main()
