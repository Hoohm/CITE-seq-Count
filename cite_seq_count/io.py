import os
import gzip
import shutil
import time
import datetime

import pandas as pd

from scipy import io
from cite_seq_count import secondsToText


def write_to_files(sparse_matrix, top_cells, ordered_tags_map, data_type, outfolder):
    """Write the umi and read sparse matrices to file in gzipped mtx format.

    Args:
        sparse_matrix (dok_matrix): Results in a sparse matrix.
        top_cells (set): Set of cells that are selected for output.
        ordered_tags_map (dict): Tags in order with indexes as values.
        data_type (string): A string definning if the data is umi or read based.
        outfolder (string): Path to the output folder.
    """
    prefix = os.path.join(outfolder, data_type + "_count")
    os.makedirs(prefix, exist_ok=True)
    io.mmwrite(os.path.join(prefix, "matrix.mtx"), sparse_matrix)
    with gzip.open(os.path.join(prefix, "barcodes.tsv.gz"), "wb") as barcode_file:
        for barcode in top_cells:
            barcode_file.write("{}\n".format(barcode).encode())
    with gzip.open(os.path.join(prefix, "features.tsv.gz"), "wb") as feature_file:
        for feature in ordered_tags_map:
            feature_file.write(
                "{}\t{}\n".format(
                    ordered_tags_map[feature]["sequence"], feature
                ).encode()
            )
    with open(os.path.join(prefix, "matrix.mtx"), "rb") as mtx_in:
        with gzip.open(os.path.join(prefix, "matrix.mtx") + ".gz", "wb") as mtx_gz:
            shutil.copyfileobj(mtx_in, mtx_gz)
    os.remove(os.path.join(prefix, "matrix.mtx"))


def write_dense(sparse_matrix, index, columns, outfolder, filename):
    """
    Writes a dense matrix in a csv format
    
    Args:
       sparse_matrix (dok_matrix): Results in a sparse matrix.
       index (list): List of TAGS
       columns (set): List of cells
       outfolder (str): Output folder
       filename (str): Filename
    """
    prefix = os.path.join(outfolder)
    os.makedirs(prefix, exist_ok=True)
    pandas_dense = pd.DataFrame(sparse_matrix.todense(), columns=columns, index=index)
    pandas_dense.to_csv(os.path.join(outfolder, filename), sep="\t")


def write_unmapped(merged_no_match, top_unknowns, outfolder, filename):
    """
    Writes a list of top unmapped sequences

    Args:
        merged_no_match (Counter): Counter of unmapped sequences
        top_unknowns (int): Number of unmapped sequences to output
        outfolder (string): Path of the output folder
        filename (string): Name of the output file
    """

    top_unmapped = merged_no_match.most_common(top_unknowns)

    with open(os.path.join(outfolder, filename), "w") as unknown_file:
        unknown_file.write("tag,count\n")
        for element in top_unmapped:
            unknown_file.write("{},{}\n".format(element[0], element[1]))


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
