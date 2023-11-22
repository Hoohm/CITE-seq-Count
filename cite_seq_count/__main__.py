#!/usr/bin/env python3.11
"""
Author: Patrick Roelli
"""
import sys
import os
import time

from cite_seq_count import preprocessing, argsparser, mapping, processing, chemistry, io


def main():
    """Main"""

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
            f"Temp folder: {args.temp_path} is not writable."
            f"Please check permissions and/or change temp folder."
        )
    if not os.access(os.path.dirname(os.path.abspath(args.outfolder)), os.W_OK):
        sys.exit(
            f"Output folder: {args.outfolder} is not writable."
            f"Please check permissions and/or change output folder."
        )

    # Get chemistry defs
    (barcode_reference, chemistry_def) = chemistry.setup_chemistry(args)

    # Load TAGs/ABs.
    parsed_tags = preprocessing.parse_tags_csv(args.tags)
    longest_tag_len = preprocessing.check_tags(parsed_tags, args.max_error)

    # Identify input file(s)
    read1_paths, read2_paths = io.get_read_paths(args.read1_path, args.read2_path)

    # Checks before chunking.
    (n_reads, r2_min_length, maximum_distance) = preprocessing.pre_run_checks(
        read1_paths=read1_paths,
        chemistry_def=chemistry_def,
        longest_tag_len=longest_tag_len,
        args=args,
    )
    (
        temp_file,
        r1_too_short,
        r2_too_short,
        total_reads,
    ) = io.write_mapping_input(
        args=args,
        read1_paths=read1_paths,
        read2_paths=read2_paths,
        r2_min_length=r2_min_length,
        chemistry_def=chemistry_def,
    )
    input_df, barcodes_df, r2_df = preprocessing.split_data_input(
        mapping_input_path=temp_file
    )
    # Remove temp file
    os.remove(temp_file)
    mapped_r2_df = mapping.map_reads_hybrid(
        r2_df=r2_df,
        parsed_tags=parsed_tags,
        maximum_distance=maximum_distance,
    )

    barcode_subset, enable_barcode_correction = preprocessing.get_barcode_subset(
        barcode_whitelist=args.filtered_barcodes,
        expected_barcodes=args.expected_barcodes,
        chemistry=chemistry_def,
        barcode_reference=barcode_reference,
        barcodes_df=barcodes_df,
    )

    # Correct cell barcodes
    if args.bc_threshold > 0 and enable_barcode_correction:
        (
            barcode_corrected_df,
            bcs_corrected,
        ) = processing.correct_barcodes_pl(
            barcodes_df=barcodes_df,
            barcode_subset_df=barcode_subset,
            hamming_distance=args.bc_threshold,
        )
    else:
        print("Skipping cell barcode correction")
        bcs_corrected = 0

    ###### HERE IT STOPS WORKING ##########

    # Create sparse matrices for reads results
    read_results_matrix = processing.generate_sparse_matrices(
        final_results=final_results,
        parsed_tags=parsed_tags,
        filtered_cells=filtered_cells,
    )
    # Write reads to file
    io.write_to_files(
        sparse_matrix=read_results_matrix,
        filtered_cells=filtered_cells,
        parsed_tags=parsed_tags,
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
            unmapped_id=len(parsed_tags),
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
            parsed_tags=parsed_tags,
            filtered_cells=clustered_cells,
        )
        # Write uncorrected cells to dense output
        io.write_dense(
            sparse_matrix=umi_clustered_matrix,
            parsed_tags=parsed_tags,
            columns=clustered_cells,
            outfolder=os.path.join(args.outfolder, "uncorrected_cells"),
            filename="dense_umis.tsv",
        )
    # Generate the UMI count matrix
    umi_results_matrix = processing.generate_sparse_matrices(
        final_results=final_results,
        parsed_tags=parsed_tags,
        filtered_cells=filtered_cells,
        umi_counts=True,
    )

    # Write umis to file
    io.write_to_files(
        sparse_matrix=umi_results_matrix,
        filtered_cells=filtered_cells,
        parsed_tags=parsed_tags,
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
            parsed_tags=parsed_tags,
            columns=filtered_cells,
            outfolder=args.outfolder,
            filename="dense_umis.tsv",
        )


if __name__ == "__main__":
    main()
