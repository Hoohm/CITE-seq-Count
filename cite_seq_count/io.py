import os
import csv
import sys
import gzip
import shutil
import time
import datetime
import tempfile

from collections import namedtuple
from itertools import islice

import pandas as pd

from scipy import io
from cite_seq_count import secondsToText


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
        if not b:
            break
        yield b


def get_n_lines(file_path):
    """
    Determines how many lines have to be processed
    depending on options and number of available lines.
    Checks that the number of lines is a multiple of 4.

    Args:
        file_path (string): Path to a fastq.gz file

    Returns:
        n_lines (int): Number of lines in the file
    """
    print("Counting number of reads in file {}".format(file_path))
    with gzip.open(file_path, "rt", encoding="utf-8", errors="ignore") as f:
        n_lines = sum(bl.count("\n") for bl in blocks(f))
    if n_lines % 4 != 0:
        sys.exit(
            "{}'s number of lines is not a multiple of 4. The file "
            "might be corrupted.\n Exiting".format(file_path)
        )
    return n_lines


def get_read_paths(read1_path, read2_path):
    """
    Splits up 2 comma-separated strings of input files into list of files
    to process. Ensures both lists are equal in length.

    Args:
        read1_path (string): Comma-separated paths to read1.fq
        read2_path (string): Comma-separated paths to read2.fq
    Returns:
        _read1_path (list(string)): list of paths to read1.fq
        _read2_path (list(string)): list of paths to read2.fq
    """
    _read1_path = read1_path.split(",")
    _read2_path = read2_path.split(",")
    if len(_read1_path) != len(_read2_path):
        sys.exit(
            "Unequal number of read1 ({}) and read2({}) files provided"
            "\n Exiting".format(len(_read1_path), len(_read2_path))
        )
    all_files = _read1_path + _read2_path
    for file_path in all_files:
        if os.path.isfile(file_path):
            if os.access(file_path, os.R_OK):
                continue
            else:
                sys.exit("{} is not readable. Exiting".format(file_path))
        else:
            sys.exit("{} does not exist. Exiting".format(file_path))

    return (_read1_path, _read2_path)


def get_csv_reader_from_path(filename, sep="\t"):
    """
    Returns a csv_reader object for a file weather it's a flat file or compressed.

    Args:
        filename: str

    Returns:
        csv_reader: The csv_reader for the file
    """
    if filename.endswith(".gz"):
        f = gzip.open(filename, mode="rt")
        csv_reader = csv.reader(f, delimiter=sep)
    else:
        f = open(filename, encoding="UTF-8")
        csv_reader = csv.reader(f, delimiter=sep)
    return csv_reader


def write_to_files(
    sparse_matrix, filtered_cells, ordered_tags, data_type, outfolder, translation_dict
):
    """Write the umi and read sparse matrices to file in gzipped mtx format.

    Args:
        sparse_matrix (dok_matrix): Results in a sparse matrix.
        filtered_cells (set): Set of cells that are selected for output.
        ordered_tags (dict): Tags in order with indexes as values.
        data_type (string): A string definning if the data is umi or read based.
        outfolder (string): Path to the output folder.
    """
    prefix = os.path.join(outfolder, data_type + "_count")
    unknown_id = 1
    os.makedirs(prefix, exist_ok=True)
    io.mmwrite(os.path.join(prefix, "matrix.mtx"), a=sparse_matrix, field="integer")
    with gzip.open(os.path.join(prefix, "barcodes.tsv.gz"), "wb") as barcode_file:
        for barcode in filtered_cells:
            if translation_dict:
                if barcode in translation_dict:
                    barcode_file.write(
                        "{}\t{}\n".format(
                            translation_dict[barcode], barcode
                        ).encode(),
                    )
                else:
                    barcode_file.write(
                        "{}\t{}\n".format(
                            "translation_not_found_{}".format(unknown_id), barcode
                        ).encode(),
                    )
                    unknown_id += 1
            else:
                barcode_file.write("{}\n".format(barcode).encode())
    with gzip.open(os.path.join(prefix, "features.tsv.gz"), "wb") as feature_file:
        for feature in ordered_tags:
            feature_file.write(
                "{}\t{}\n".format(feature.sequence, feature.name).encode()
            )
        if data_type == "read":
            feature_file.write("{}\t{}\n".format("UNKNOWN", "unmapped").encode())
    with open(os.path.join(prefix, "matrix.mtx"), "rb") as mtx_in:
        with gzip.open(os.path.join(prefix, "matrix.mtx") + ".gz", "wb") as mtx_gz:
            shutil.copyfileobj(mtx_in, mtx_gz)
    os.remove(os.path.join(prefix, "matrix.mtx"))


def write_dense(sparse_matrix, ordered_tags, columns, outfolder, filename):
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
    index = []
    for tag in ordered_tags:
        index.append(tag.name)
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
    ordered_tags,
    umis_corrected,
    bcs_corrected,
    bad_cells,
    R1_too_short,
    R2_too_short,
    args,
    chemistry_def,
    maximum_distance,
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
    total_too_short = R1_too_short + R2_too_short
    total_mapped = total_reads - total_unmapped - total_too_short

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
  R1_too_short: {}
  R2_too_short: {}
Uncorrected cells: {}
Correction:
  Cell barcodes collapsing threshold: {}
  Cell barcodes corrected: {}
  UMI collapsing threshold: {}
  UMIs corrected: {}
Run parameters:
  Read1_paths: {}
  Read2_paths: {}
  Cell barcode:
    First position: {}
    Last position: {}
  UMI barcode:
    First position: {}
    Last position: {}
  Expected cells: {}
  Tags max errors: {}
  Start trim: {}
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
                chemistry_def.cell_barcode_start,
                chemistry_def.cell_barcode_end,
                chemistry_def.umi_barcode_start,
                chemistry_def.umi_barcode_end,
                args.expected_cells,
                maximum_distance,
                chemistry_def.R2_trim_start,
            )
        )


def write_chunks_to_disk(
    args,
    read1_paths,
    read2_paths,
    R2_min_length,
    n_reads_per_chunk,
    chemistry_def,
    ordered_tags,
    maximum_distance,
):
    """
    Writes chunked files of reads to disk and prepares parallel
    processing queue parameters.

    Args:
        args(argparse): All parsed arguments.
        read1_paths (list): List of R1 fastq.gz paths.
        read2_paths (list): List of R2 fastq.gz paths.
        R2_min_length (int):  Minimum length of read2 sequences.
        n_reads_per_chunk (int): How many reads per chunk.
        chemistry_def (namedtuple): Hols all the information about the chemistry definition.
        ordered_tags (list): List of namedtuple tags.
        maximum_distance (int): Maximum hamming distance for mapping.
    """
    mapping_input = namedtuple(
        "mapping_input",
        ["filename", "tags", "debug", "maximum_distance", "sliding_window"],
    )

    print("Writing chunks to disk")

    num_chunk = 0
    if not args.chunk_size:
        chunk_size = round(n_reads_per_chunk / args.n_threads)
    else:
        chunk_size = args.chunk_size
    temp_path = os.path.abspath(args.temp_path)
    input_queue = []
    temp_files = []
    R1_too_short = 0
    R2_too_short = 0
    total_reads = 0
    total_reads_written = 0
    enough_reads = False

    barcode_slice = slice(
        chemistry_def.cell_barcode_start - 1, chemistry_def.cell_barcode_end
    )
    umi_slice = slice(
        chemistry_def.umi_barcode_start - 1, chemistry_def.umi_barcode_end
    )

    chunked_file_object = tempfile.NamedTemporaryFile(
        "w", dir=temp_path, suffix="_csc", delete=False
    )
    # chunked_file_object = open(temp_file, "w")
    temp_files.append(chunked_file_object.name)
    reads_written = 0

    for read1_path, read2_path in zip(read1_paths, read2_paths):

        if enough_reads:
            break
        print("Reading reads from files: {}, {}".format(read1_path, read2_path))
        with gzip.open(read1_path, "rt") as textfile1, gzip.open(
            read2_path, "rt"
        ) as textfile2:
            secondlines = islice(zip(textfile1, textfile2), 1, None, 4)

            for read1, read2 in secondlines:
                total_reads += 1

                read1 = read1.strip()
                if len(read1) < chemistry_def.umi_barcode_end:
                    R1_too_short += 1
                    # The entire read is skipped
                    continue
                if len(read2) < R2_min_length:
                    R2_too_short += 1
                    # The entire read is skipped
                    continue

                read1_sliced = read1[
                    chemistry_def.cell_barcode_start - 1 : chemistry_def.umi_barcode_end
                ]

                read2_sliced = read2[
                    chemistry_def.R2_trim_start : (
                        R2_min_length + chemistry_def.R2_trim_start
                    )
                ]
                chunked_file_object.write(
                    "{},{},{}\n".format(
                        read1_sliced[barcode_slice],
                        read1_sliced[umi_slice],
                        read2_sliced,
                    )
                )

                reads_written += 1
                total_reads_written += 1
                if reads_written % chunk_size == 0 and reads_written != 0:
                    # We have enough reads in this chunk, open a new one
                    chunked_file_object.close()
                    input_queue.append(
                        mapping_input(
                            filename=chunked_file_object.name,
                            tags=ordered_tags,
                            debug=args.debug,
                            maximum_distance=maximum_distance,
                            sliding_window=args.sliding_window,
                        )
                    )
                    if total_reads_written == n_reads_per_chunk:
                        enough_reads = True
                        chunked_file_object.close()
                        break
                    num_chunk += 1
                    chunked_file_object = tempfile.NamedTemporaryFile(
                        "w", dir=temp_path, suffix="_csc", delete=False
                    )
                    temp_files.append(chunked_file_object.name)
                    reads_written = 0
                if total_reads_written == n_reads_per_chunk:
                    enough_reads = True
                    chunked_file_object.close()
                    break
    if not enough_reads:
        chunked_file_object.close()
        input_queue.append(
            mapping_input(
                filename=chunked_file_object.name,
                tags=ordered_tags,
                debug=args.debug,
                maximum_distance=maximum_distance,
                sliding_window=args.sliding_window,
            )
        )
    return (
        input_queue,
        temp_files,
        R1_too_short,
        R2_too_short,
        total_reads,
    )
