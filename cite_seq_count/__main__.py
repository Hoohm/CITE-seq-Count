#!/usr/bin/env python3.11
"""
Author: Patrick Roelli
"""
import sys
import os
import time

from scipy import constants

from cite_seq_count import (
    preprocessing,
    argsparser,
    mapping,
    processing,
    chemistry,
    io,
    constants,
)


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
            f"Temp folder: {args.temp_path} is not writeable."
            f"Please check permissions and/or change temp folder."
        )
    if not os.access(os.path.dirname(os.path.abspath(args.outfolder)), os.W_OK):
        sys.exit(
            f"Output folder: {args.outfolder} is not writeable."
            f"Please check permissions and/or change output folder."
        )

    # Get chemistry defs
    barcode_reference, chemistry_def = chemistry.setup_chemistry(args)
    if args.subset_path is not None:
        barcode_subset = preprocessing.parse_barcode_file(
            filename=args.subset_path,
            barcode_length=chemistry_def.barcode_length,
            required_header=[constants.REFERENCE_COLUMN],
        )
    else:
        barcode_subset = None

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
        arguments=args,
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
    main_df, barcodes_df, r2_df = preprocessing.split_data_input(
        mapping_input_path=temp_file
    )
    # Remove temp file
    os.remove(temp_file)
    mapped_r2_df, unmapped_r2_df = mapping.map_reads_hybrid(
        r2_df=r2_df,
        parsed_tags=parsed_tags,
        maximum_distance=maximum_distance,
    )
    unmapped_df = processing.summarise_unmapped_df(
        main_df=main_df, unmapped_r2_df=unmapped_r2_df
    )
    barcode_subset, enable_barcode_correction = preprocessing.get_barcode_subset(
        barcode_subset=barcode_subset,
        n_barcodes=args.expected_barcodes,
        chemistry=chemistry_def,
        barcode_reference=barcode_reference,
        barcodes_df=barcodes_df,
    )

    # Correct cell barcodes
    if args.bc_threshold > 0 and enable_barcode_correction:
        (
            barcodes_df,
            n_bcs_corrected,
            mapped_barcodes,
        ) = processing.correct_barcodes_pl(
            barcodes_df=barcodes_df,
            barcode_subset_df=barcode_subset,
            hamming_distance=args.bc_threshold,
        )
        main_df = processing.update_main_df(
            main_df=main_df, mapped_barcodes=mapped_barcodes
        )
    else:
        print("Skipping cell barcode correction")
        n_bcs_corrected = 0
    read_counts = processing.generate_mtx_counts(
        main_df=main_df,
        barcode_subset=barcodes_df,
        mapped_r2_df=mapped_r2_df,
        data_type="read",
    )
    # Write reads to file
    io.write_data_to_mtx(
        main_df=read_counts,
        tags_df=parsed_tags,
        subset_df=barcode_subset,
        data_type="read",
        outpath=args.outfolder,
    )
    # TODO: add clustered cells filter: Max UMIs per cell per feature: 20000
    print("UMI correction not implemented yet")
    # Don't correct
    umis_corrected = 0
    clustered_cells = []
    # TODO: Write out to mtx and csv clustered cells

    # Generate the UMI count matrix
    umi_counts = processing.generate_mtx_counts(
        main_df=main_df,
        barcode_subset=barcodes_df,
        mapped_r2_df=mapped_r2_df,
        data_type="umi",
    )

    # Write umis to file
    io.write_data_to_mtx(
        main_df=read_counts,
        tags_df=parsed_tags,
        subset_df=barcode_subset,
        data_type="read",
        outpath=args.outfolder,
    )
    io.write_data_to_mtx(
        main_df=umi_counts,
        tags_df=parsed_tags,
        subset_df=barcode_subset,
        data_type="umi",
        outpath=args.outfolder,
    )

    # TODO: Write unmapped sequences
    # TODO: rewrite reporting
    # Create report and write it to disk
    io.create_report(
        total_reads=total_reads,
        unmapped=unmapped_df,
        version=argsparser.get_package_version(),
        start_time=start_time,
        umis_corrected=umis_corrected,
        bcs_corrected=n_bcs_corrected,
        bad_cells=clustered_cells,
        r1_too_short=r1_too_short,
        r2_too_short=r2_too_short,
        args=args,
        chemistry_def=chemistry_def,
        maximum_distance=maximum_distance,
    )


if __name__ == "__main__":
    main()
