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
import gzip
import requests

from itertools import islice
from argparse import ArgumentParser, ArgumentTypeError, RawTextHelpFormatter

from collections import OrderedDict, Counter, defaultdict, namedtuple
from multiprocess import cpu_count, Pool, Queue, JoinableQueue, Process
from tempfile import NamedTemporaryFile

from cite_seq_count import preprocessing
from cite_seq_count import processing
from cite_seq_count import database
from cite_seq_count import io
from cite_seq_count import secondsToText

version = pkg_resources.require("cite_seq_count")[0].version


def chunk_size_limit(arg):
    """Validates chunk_size limits"""
    max_size = 2147483647
    try:
        f = int(arg)
    except ValueError:
        raise ArgumentTypeError("Must be am int")
    if f < 1 or f > max_size:
        raise ArgumentTypeError(
            "Argument must be < " + str(max_size) + "and > " + str(1)
        )
    return f


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
    barcodes.add_argument("--chemistry", type=str, required=False, default=False)
    if "--chemistry" not in sys.argv:
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
    # Cells group
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
    if "--chemistry" not in sys.argv:
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

        cells.add_argument(
            "--translation",
            required=False,
            type=str,
            help="A csv file containing the mapping between two sets of cell barcode list.\n"
            "A required header such as the reference is named whitelist. Example:\n\n"
            "\twhitelist,trasnlation\n"
            "\tAAACCCAAGAAACACT,AAACCCATCAAACACT\n"
            "\tAAACCCAAGAAACCAT,AAACCCATCAAACCAT\n"
            "\nThe output matrix will possess both cell barcode IDs",
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

    # Parallel group.
    parallel = parser.add_argument_group(
        "Parallelization options",
        description=("Options for performance on parallelization"),
    )
    # Remaining arguments.
    parallel.add_argument(
        "-T",
        "--threads",
        required=False,
        type=int,
        dest="n_threads",
        default=cpu_count(),
        help="How many threads are to be used for running the program",
    )
    parallel.add_argument(
        "-C",
        "--chunk_size",
        required=False,
        type=chunk_size_limit,
        dest="chunk_size",
        help="How many reads should be sent to a child process at a time",
    )
    parallel.add_argument(
        "--temp_path",
        required=False,
        type=str,
        dest="temp_path",
        default="",
        help="Temp folder for chunk creation specification. Useful when using a cluster with a scratch folder",
    )

    # Global group
    parser.add_argument(
        "-n",
        "--first_n",
        required=False,
        type=int,
        dest="first_n",
        default=float("inf"),
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
    total_reads,
    reads_per_cell,
    no_match,
    version,
    start_time,
    ordered_tags_map,
    umis_corrected,
    bcs_corrected,
    bad_cells,
    R1_too_short,
    R2_too_short,
    args,
    chemistry_def,
):
    """
    Creates a report with details about the run in a yaml format.
    Args:
        total_reads (int): Number of reads that have been processed.
        reads_matrix (scipy.sparse.dok_matrix): A sparse matrix continining read counts.
        no_match (Counter): Counter of unmapped tags.
        version (string): CITE-seq-Count package version.
        start_time (time): Start time of the run.
        args (arg_parse): Arguments provided by the user.

    """
    total_unmapped = sum(no_match.values())
    total_mapped = total_reads - total_unmapped
    total_too_short = total_reads - total_unmapped - total_mapped
    too_short_perc = round((total_too_short / total_reads) * 100)
    mapped_perc = round((total_mapped / total_reads) * 100)
    unmapped_perc = round((total_unmapped / total_reads) * 100)

    with open(os.path.join(args.outfolder, "run_report.yaml"), "w") as report_file:
        report_file.write(
            """Date: {}
Running time: {}
CITE-seq-Count Version: {}
Reads processed: {}
Percentage mapped: {}
Percentage unmapped: {}
Percentage too short: {}
\tR1_too_short: {}
\tR2_too_short: {}
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
                int(total_reads),
                mapped_perc,
                unmapped_perc,
                too_short_perc,
                R1_too_short,
                R2_too_short,
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
                chemistry_def.umi_barcode_end,
                args.expected_cells,
                args.max_error,
                chemistry_def.R2_trim_start,
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
    temp_path = os.path.abspath(args.temp_path)
    assert os.access(temp_path, os.W_OK)

    # Get chemistry defs
    if args.chemistry:
        chemistry_def = database.get_chemistry_definition(args.chemistry)
        with requests.get(chemistry_def.whitelist_path) as r:
            r.raise_for_status()
            ###### TODO deal with the download of the whitelist from remote
            with NamedTemporaryFile() as temp_local_whitelist:
                temp_local_whitelist.write(r.content)
                (whitelist, args.bc_threshold) = preprocessing.parse_whitelist_csv(
                    filename=temp_local_whitelist,
                    barcode_length=chemistry_def.cell_barcode_end
                    - chemistry_def.cell_barcode_start
                    + 1,
                    collapsing_threshold=args.bc_threshold,
                )
    else:
        chemistry_def = database.create_chemistry_definition(args)
        if args.whitelist:
            print("Loading whitelist")
            (whitelist, args.bc_threshold) = preprocessing.parse_whitelist_csv(
                filename=args.whitelist,
                barcode_length=args.cb_last - args.cb_first + 1,
                collapsing_threshold=args.bc_threshold,
            )
        else:
            whitelist = False

    # Load TAGs/ABs.
    ab_map = preprocessing.parse_tags_csv(args.tags)
    ordered_tags_map, longest_tag_len = preprocessing.check_tags(ab_map, args.max_error)
    named_tuples_tags_map = preprocessing.convert_to_named_tuple(
        ordered_tags=ordered_tags_map
    )
    # Identify input file(s)
    read1_paths, read2_paths = preprocessing.get_read_paths(
        args.read1_path, args.read2_path
    )

    # preprocessing and processing occur in separate loops so the program can crash earlier if
    # one of the inputs is not valid.
    read1_lengths = []
    read2_lengths = []
    total_reads = 0

    for read1_path, read2_path in zip(read1_paths, read2_paths):
        n_lines = preprocessing.get_n_lines(read1_path)
        total_reads += n_lines / 4
        # Get reads length. So far, there is no validation for Read2.
        read1_lengths.append(preprocessing.get_read_length(read1_path))
        read2_lengths.append(preprocessing.get_read_length(read2_path))
        # Check Read1 length against CELL and UMI barcodes length.
        (barcode_slice, umi_slice, _) = preprocessing.check_barcodes_lengths(
            read1_lengths[-1],
            chemistry_def.cell_barcode_start,
            chemistry_def.cell_barcode_end,
            chemistry_def.umi_barcode_start,
            chemistry_def.umi_barcode_end,
        )

    # Ensure all files have the same input length
    # if len(set(read1_lengths)) != 1:
    # sys.exit('Input barcode fastqs (read1) do not all have same length.\nExiting')
    # if len(set(read2_lengths)) != 1:
    # sys.exit('Input barcode fastqs (read2) do not all have same length.\nExiting')

    # Define R2_lenght to reduce amount of data to transfer to childrens
    if args.sliding_window:
        R2_max_length = read2_lengths[0]
    else:
        R2_max_length = longest_tag_len
    # Initialize the counts dicts that will be generated from each input fastq pair
    final_results = defaultdict(lambda: defaultdict(Counter))
    umis_per_cell = Counter()
    reads_per_cell = Counter()
    merged_no_match = Counter()
    number_of_samples = len(read1_paths)

    # Print a statement if multiple files are run.
    if number_of_samples != 1:
        print("Detected {} files to run on.".format(number_of_samples))
    input_queue = []
    mapping_input = namedtuple(
        "mapping_input",
        ["filename", "tags", "debug", "maximum_distance", "sliding_window"],
    )

    print("Writing chunks to disk")
    reads_count = 0
    num_chunks = 0
    if args.chunk_size:
        chunk_size = args.chunk_size
    else:
        chunk_size = round(total_reads / args.n_threads) + 1
    temp_files = []
    R1_too_short = 0
    R2_too_short = 0
    for read1_path, read2_path in zip(read1_paths, read2_paths):
        print("Reading reads from files: {}, {}".format(read1_path, read2_path))
        with gzip.open(read1_path, "rt") as textfile1, gzip.open(
            read2_path, "rt"
        ) as textfile2:
            secondlines = islice(zip(textfile1, textfile2), 1, None, 4)
            temp_filename = os.path.join(temp_path, "temp_{}".format(num_chunks))
            chunked_file_object = open(temp_filename, "w")
            temp_files.append(os.path.abspath(temp_filename))
            for read1, read2 in secondlines:

                read1 = read1.strip()
                if len(read1) < chemistry_def.umi_barcode_end:
                    R1_too_short += 1
                    # The entire read is skipped
                    continue
                read1_sliced = read1[0 : chemistry_def.umi_barcode_end]
                if len(read2) < R2_max_length:
                    R2_too_short += 1
                    # The entire read is skipped
                    continue

                read2_sliced = read2[
                    chemistry_def.R2_trim_start : (
                        R2_max_length + chemistry_def.R2_trim_start
                    )
                ]
                chunked_file_object.write(
                    "{},{},{}\n".format(
                        read1_sliced[barcode_slice],
                        read1_sliced[umi_slice],
                        read2_sliced,
                    )
                )
                reads_count += 1
                if reads_count % chunk_size == 0:
                    input_queue.append(
                        mapping_input(
                            filename=temp_filename,
                            tags=named_tuples_tags_map,
                            debug=args.debug,
                            maximum_distance=args.max_error,
                            sliding_window=args.sliding_window,
                        )
                    )
                    num_chunks += 1
                    chunked_file_object.close()
                    temp_filename = "temp_{}".format(num_chunks)
                    chunked_file_object = open(temp_filename, "w")
                    temp_files.append(os.path.abspath(temp_filename))
                if reads_count >= args.first_n:
                    total_reads = args.first_n
                    break

            input_queue.append(
                mapping_input(
                    filename=temp_filename,
                    tags=named_tuples_tags_map,
                    debug=args.debug,
                    maximum_distance=args.max_error,
                    sliding_window=args.sliding_window,
                )
            )
            chunked_file_object.close()

    print("Started mapping")
    parallel_results = []
    pool = Pool(processes=args.n_threads)
    errors = []
    mapping = pool.map_async(
        processing.map_reads,
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
    (
        final_results,
        umis_per_cell,
        reads_per_cell,
        merged_no_match,
    ) = processing.merge_results(parallel_results=parallel_results[0])

    del parallel_results

    # Check if 99% of the reads are unmapped.
    processing.check_unmapped(
        no_match=merged_no_match,
        total_reads=total_reads,
        start_trim=chemistry_def.R2_trim_start,
    )
    # Delete temp_files
    for file_path in temp_files:
        if os.path.exists(file_path):
            os.remove(file_path)
        else:
            print("Could not find file: {}".format(file_path))

    # Correct cell barcodes
    if args.bc_threshold != 0:
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
            if not whitelist:
                (
                    final_results,
                    umis_per_cell,
                    bcs_corrected,
                ) = processing.correct_cells(
                    final_results=final_results,
                    reads_per_cell=reads_per_cell,
                    umis_per_cell=umis_per_cell,
                    expected_cells=args.expected_cells,
                    collapsing_threshold=args.bc_threshold,
                    ab_map=named_tuples_tags_map,
                )
            else:
                (
                    final_results,
                    umis_per_cell,
                    bcs_corrected,
                ) = processing.correct_cells_whitelist(
                    final_results=final_results,
                    umis_per_cell=umis_per_cell,
                    whitelist=whitelist,
                    collapsing_threshold=args.bc_threshold,
                    ab_map=named_tuples_tags_map,
                )
    else:
        print("Skipping cell barcode correction")
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
                for TAG in named_tuples_tags_map:
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
        input_queue = []

        umi_correction_input = namedtuple(
            "umi_correction_input", ["cells", "collapsing_threshold", "max_umis"]
        )
        cells = {}
        n_cells = 0
        num_chunks = 0

        cell_batch_size = round(len(top_cells) / args.n_threads) + 1
        for cell in top_cells:
            cells[cell] = final_results[cell]
            n_cells += 1
            if n_cells % cell_batch_size == 0:
                input_queue.append(
                    umi_correction_input(
                        cells=cells,
                        collapsing_threshold=args.umi_threshold,
                        max_umis=20000,
                    )
                )
                cells = {}
                num_chunks += 1
        input_queue.append(
            umi_correction_input(
                cells=cells, collapsing_threshold=args.umi_threshold, max_umis=20000
            )
        )

        pool = Pool(processes=args.n_threads)
        errors = []
        parallel_results = []
        correct_umis = pool.map_async(
            processing.correct_umis,
            input_queue,
            callback=parallel_results.append,
            error_callback=errors.append,
        )

        correct_umis.wait()
        pool.close()
        pool.join()

        if len(errors) != 0:
            for error in errors:
                print(error)

        final_results = {}
        umis_corrected = 0
        aberrant_cells = set()

        for chunk in parallel_results[0]:
            (temp_results, temp_umis, temp_aberrant_cells) = chunk
            final_results.update(temp_results)
            umis_corrected += temp_umis
            aberrant_cells.update(temp_aberrant_cells)

    if len(aberrant_cells) > 0:
        # Remove aberrant cells from the top cells
        for cell_barcode in aberrant_cells:
            top_cells.remove(cell_barcode)

        # Create sparse aberrant cells matrix
        (umi_aberrant_matrix, _) = processing.generate_sparse_matrices(
            final_results=final_results,
            ordered_tags_map=ordered_tags_map,
            top_cells=aberrant_cells,
        )

        # Write uncorrected cells to dense output
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
        total_reads=total_reads,
        reads_per_cell=reads_per_cell,
        no_match=merged_no_match,
        version=version,
        start_time=start_time,
        ordered_tags_map=ordered_tags_map,
        umis_corrected=umis_corrected,
        bcs_corrected=bcs_corrected,
        bad_cells=aberrant_cells,
        R1_too_short=R1_too_short,
        R2_too_short=R2_too_short,
        args=args,
        chemistry_def=chemistry_def,
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


if __name__ == "__main__":
    main()
