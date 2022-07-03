#!/usr/bin/env python3.6
"""
Author: Patrick Roelli
"""
import sys
import os
import time

from cite_seq_count import preprocessing, argsparser, mapping, processing, chemistry, io


def main():
    """Main function"""

    start_time = time.time()
    parser = argsparser.get_args()
    if not sys.argv[1:]:
        parser.print_help(file=sys.stderr)
        sys.exit(2)

    # Parse arguments.
    args = parser.parse_args()
    # Check a few path before doing anything
    if not os.access(args.temp_path, os.W_OK):
        sys.exit(
            "Temp folder: {} is not writable. Please check permissions and/or change temp folder.".format(
                args.temp_path
            )
        )
    if not os.access(os.path.dirname(os.path.abspath(args.outfolder)), os.W_OK):
        sys.exit(
            "Output folder: {} is not writable. Please check permissions and/or change output folder.".format(
                args.outfolder
            )
        )

    # Get chemistry defs
    (translation_dict, chemistry_def) = chemistry.setup_chemistry(args)

    # Load TAGs/ABs.
    ab_map = preprocessing.parse_tags_csv(args.tags)
    ordered_tags, longest_tag_len = preprocessing.check_tags(ab_map, args.max_error)

    # Identify input file(s)
    read1_paths, read2_paths = io.get_read_paths(args.read1_path, args.read2_path)

    # Check filtered input list
    # If a translation is given, will return the translated version
    filtered_cells = preprocessing.get_filtered_list(
        args=args, chemistry=chemistry_def, translation_dict=translation_dict
    )
    # Checks before chunking.
    (n_reads, r2_min_length, maximum_distance) = preprocessing.pre_run_checks(
        read1_paths=read1_paths,
        chemistry_def=chemistry_def,
        longest_tag_len=longest_tag_len,
        args=args,
    )

    # Chunk the data to disk before mapping
    (
        input_queue,
        temp_files,
        r1_too_short,
        r2_too_short,
        total_reads,
    ) = io.write_chunks_to_disk(
        args=args,
        read1_paths=read1_paths,
        read2_paths=read2_paths,
        r2_min_length=r2_min_length,
        n_reads_per_chunk=n_reads,
        chemistry_def=chemistry_def,
        ordered_tags=ordered_tags,
        maximum_distance=maximum_distance,
    )
    # Map the data
    (final_results, umis_per_cell, reads_per_cell, merged_no_match) = mapping.map_data(
        input_queue=input_queue, unmapped_id=len(ordered_tags), args=args
    )

    # Check if 99% of the reads are unmapped.
    mapping.check_unmapped(
        no_match=merged_no_match,
        too_short=r1_too_short + r2_too_short,
        total_reads=total_reads,
        start_trim=chemistry_def.r2_trim_start,
    )

    # Remove temp chunks
    for file_path in temp_files:
        os.remove(file_path)
    # Check that we have a filtered cell list to work on
    filtered_cells = processing.check_filtered_cells(
        filtered_cells=filtered_cells,
        expected_cells=args.expected_cells,
        umis_per_cell=umis_per_cell,
    )

    # Correct cell barcodes
    if args.bc_threshold > 0:
        (
            final_results,
            umis_per_cell,
            bcs_corrected,
        ) = processing.run_cell_barcode_correction(
            final_results=final_results,
            umis_per_cell=umis_per_cell,
            ordered_tags=ordered_tags,
            filtered_set=filtered_cells,
            args=args,
        )
    else:
        print("Skipping cell barcode correction")
        bcs_corrected = 0

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
        translation_dict=translation_dict,
    )

    # UMI correction
    if args.umi_threshold != 0:
        # Correct UMIS
        (
            final_results,
            umis_corrected,
            clustered_cells,
        ) = processing.run_umi_correction(
            final_results=final_results,
            filtered_cells=filtered_cells,
            unmapped_id=len(ordered_tags),
            args=args,
        )
    else:
        # Don't correct
        umis_corrected = 0
        clustered_cells = []

    if len(clustered_cells) > 0:
        # Remove clustered cells from the top cells
        for cell_barcode in clustered_cells:
            filtered_cells.remove(cell_barcode)

        # Create sparse clustered cells matrix
        umi_clustered_matrix = processing.generate_sparse_matrices(
            final_results=final_results,
            ordered_tags=ordered_tags,
            filtered_cells=clustered_cells,
        )
        # Write uncorrected cells to dense output
        io.write_dense(
            sparse_matrix=umi_clustered_matrix,
            ordered_tags=ordered_tags,
            columns=clustered_cells,
            outfolder=os.path.join(args.outfolder, "uncorrected_cells"),
            filename="dense_umis.tsv",
        )
    # Generate the UMI count matrix
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
        translation_dict=translation_dict,
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
        no_match=merged_no_match,
        version=argsparser.get_package_version(),
        start_time=start_time,
        umis_corrected=umis_corrected,
        bcs_corrected=bcs_corrected,
        bad_cells=clustered_cells,
        r1_too_short=r1_too_short,
        r2_too_short=r2_too_short,
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
