#!/usr/bin/env python3.6
"""
Author: Patrick Roelli
"""
import sys
import os
import logging
import gzip
import requests
import time

from collections import OrderedDict, Counter, defaultdict, namedtuple

# pylint: disable=no-name-in-module
from multiprocess import Pool, Queue, JoinableQueue, Process

from cite_seq_count import preprocessing
from cite_seq_count import processing
from cite_seq_count import chemistry
from cite_seq_count import io
from cite_seq_count import secondsToText
from cite_seq_count import argsparser


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
    parser = argsparser.get_args()
    if not sys.argv[1:]:
        parser.print_help(file=sys.stderr)
        sys.exit(2)

    # Parse arguments.
    args = parser.parse_args()
    assert os.access(args.temp_path, os.W_OK)

    # Get chemistry defs
    (whitelist, chemistry_def) = chemistry.setup_chemistry(args)

    # Load TAGs/ABs.
    ab_map = preprocessing.parse_tags_csv(args.tags)
    ordered_tags, longest_tag_len = preprocessing.check_tags(ab_map, args.max_error)
    # ordered_tags = preprocessing.convert_to_named_tuple(ordered_tags=ordered_tags)
    # Identify input file(s)
    read1_paths, read2_paths = preprocessing.get_read_paths(
        args.read1_path, args.read2_path
    )

    # preprocessing and processing occur in separate loops so the program can crash earlier if
    # one of the inputs is not valid.
    read1_lengths = []
    read2_lengths = []
    total_reads = 0

    for read1_path in read1_paths:
        n_lines = preprocessing.get_n_lines(read1_path)
        total_reads += n_lines / 4
        # Get reads length. So far, there is no validation for Read2.
        read1_lengths.append(preprocessing.get_read_length(read1_path))
        # read2_lengths.append(preprocessing.get_read_length(read2_path))
        # Check Read1 length against CELL and UMI barcodes length.
        preprocessing.check_barcodes_lengths(
            read1_lengths[-1],
            chemistry_def.cell_barcode_start,
            chemistry_def.cell_barcode_end,
            chemistry_def.umi_barcode_start,
            chemistry_def.umi_barcode_end,
        )

    # Get all reads or only top N?
    if args.first_n < float("inf"):
        n_reads = args.first_n
    else:
        n_reads = total_reads

    # Define R2_lenght to reduce amount of data to transfer to childrens
    number_of_samples = len(read1_paths)

    # Print a statement if multiple files are run.
    if number_of_samples != 1:
        print("Detected {} pairs of files to run on.".format(number_of_samples))

    if args.sliding_window:
        R2_min_length = read2_lengths[0]
        maximum_distance = 0
    else:
        R2_min_length = longest_tag_len
        maximum_distance = args.max_error

    (
        input_queue,
        temp_files,
        R1_too_short,
        R2_too_short,
        total_reads,
    ) = io.write_chunks_to_disk(
        args=args,
        read1_paths=read1_paths,
        read2_paths=read2_paths,
        R2_min_length=R2_min_length,
        n_reads_to_chunk=n_reads,
        chemistry_def=chemistry_def,
        ordered_tags=ordered_tags,
        maximum_distance=maximum_distance,
    )
    # Initialize the counts dicts that will be generated from each input fastq pair
    final_results = defaultdict(lambda: defaultdict(Counter))
    umis_per_cell = Counter()
    reads_per_cell = Counter()
    merged_no_match = Counter()

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
        too_short=R1_too_short + R2_too_short,
        total_reads=total_reads,
        start_trim=chemistry_def.R2_trim_start,
    )
    # Delete temp_files
    # exit()
    for file_path in temp_files:
        os.remove(file_path)

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
                    ab_map=ordered_tags,
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
                    ab_map=ordered_tags,
                )
    else:
        print("Skipping cell barcode correction")
        bcs_corrected = 0

    # If given, use whitelist for top cells
    top_cells_tuple = umis_per_cell.most_common(args.expected_cells * 10)
    if whitelist:
        # Add potential missing cell barcodes.
        # for missing_cell in whitelist:
        #     if missing_cell in final_results:
        #         continue
        #     else:
        #         final_results[missing_cell] = dict()
        #         for TAG in ordered_tags:
        #             final_results[missing_cell][TAG.safe_name] = Counter()
        #         filtered_cells.add(missing_cell)
        top_cells = set([pair[0] for pair in top_cells_tuple])
        filtered_cells = set()
        for cell in top_cells:
            if cell in whitelist:
                filtered_cells.add(cell)
    else:
        # Select top cells based on total umis per cell
        filtered_cells = set([pair[0] for pair in top_cells_tuple])

    # Create sparse matrices for reads results
    read_results_matrix = processing.generate_sparse_matrices(
        final_results=final_results,
        ordered_tags=ordered_tags,
        filtered_cells=filtered_cells,
    )
    # Write reads to file
    io.write_to_files(
        sparse_matrix=read_results_matrix,
        filtered_cells=filtered_cells,
        ordered_tags=ordered_tags,
        data_type="read",
        outfolder=args.outfolder,
    )

    # UMI correction
    if args.umi_threshold != 0:
        # Correct UMIS
        input_queue = []

        umi_correction_input = namedtuple(
            "umi_correction_input", ["cells", "collapsing_threshold", "max_umis"]
        )
        cells = {}
        n_cells = 0
        num_chunks = 0
        print("preparing UMI correction jobs")
        cell_batch_size = round(len(filtered_cells) / args.n_threads) + 1
        for cell in filtered_cells:
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
    else:
        # Don't correct
        umis_corrected = 0
        aberrant_cells = set()

    if len(aberrant_cells) > 0:
        # Remove aberrant cells from the top cells
        for cell_barcode in aberrant_cells:
            filtered_cells.remove(cell_barcode)

        # Create sparse aberrant cells matrix
        umi_aberrant_matrix = processing.generate_sparse_matrices(
            final_results=final_results,
            ordered_tags=ordered_tags,
            filtered_cells=aberrant_cells,
        )

        # Write uncorrected cells to dense output
        io.write_dense(
            sparse_matrix=umi_aberrant_matrix,
            ordered_tags=ordered_tags,
            columns=aberrant_cells,
            outfolder=os.path.join(args.outfolder, "uncorrected_cells"),
            filename="dense_umis.tsv",
        )
    # delete the last element (unmapped)
    ordered_tags.pop()
    umi_results_matrix = processing.generate_sparse_matrices(
        final_results=final_results,
        ordered_tags=ordered_tags,
        filtered_cells=filtered_cells,
        umi_counts=True,
    )

    # Write umis to file
    io.write_to_files(
        sparse_matrix=umi_results_matrix,
        filtered_cells=filtered_cells,
        ordered_tags=ordered_tags,
        data_type="umi",
        outfolder=args.outfolder,
    )

    # Write unmapped sequences
    if len(merged_no_match) > 0:
        io.write_unmapped(
            merged_no_match=merged_no_match,
            top_unknowns=args.unknowns_top,
            outfolder=args.outfolder,
            filename=args.unmapped_file,
        )

    # Create report and write it to disk
    io.create_report(
        total_reads=total_reads,
        reads_per_cell=reads_per_cell,
        no_match=merged_no_match,
        version=argsparser.get_package_version(),
        start_time=start_time,
        ordered_tags=ordered_tags,
        umis_corrected=umis_corrected,
        bcs_corrected=bcs_corrected,
        bad_cells=aberrant_cells,
        R1_too_short=R1_too_short,
        R2_too_short=R2_too_short,
        args=args,
        chemistry_def=chemistry_def,
        maximum_distance=maximum_distance,
    )

    # Write dense matrix to disk if requested
    if args.dense:
        print("Writing dense format output")
        io.write_dense(
            sparse_matrix=umi_results_matrix,
            ordered_tags=ordered_tags,
            columns=filtered_cells,
            outfolder=args.outfolder,
            filename="dense_umis.tsv",
        )


if __name__ == "__main__":
    main()
