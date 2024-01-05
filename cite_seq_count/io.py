"""Handle io operations"""
import os
import csv
import sys
import gzip
import shutil
import time
import datetime
import tempfile
import json

from argparse import Namespace
from itertools import islice
from collections import Counter
from typing import Tuple, TextIO
from pathlib import Path
from os import access, R_OK

import scipy
import pkg_resources
import yaml
import pandas as pd
import polars as pl
from scipy import io
from cite_seq_count import secondsToText
from cite_seq_count.constants import (
    FEATURE_NAME_COLUMN,
    BARCODE_COLUMN,
    BARCODE_ID_COLUMN,
    FEATURE_ID_COLUMN,
    COUNT_COLUMN,
    MTX_HEADER,
    FEATURES_MTX,
    BARCODE_MTX,
    MATRIX_MTX,
    TEMP_MTX,
    UNMAPPED_NAME,
    SEQUENCE_COLUMN,
    SUBSET_COLUMN,
)

JSON_REPORT_PATH = pkg_resources.resource_filename(__name__, "templates/report.json")


def compress_file(in_file: Path, out_file: Path):
    with open(in_file, "rb") as f_in:
        with gzip.open(out_file, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
    os.remove(in_file)


def blocks(file: TextIO, size: int = 65536):
    """
    A fast way of counting the lines of a large file.
    Ref:
        https://stackoverflow.com/a/9631635/9178565

    Args:
        file (io.handler): A file handler
        size (int): Block size
    Returns:
        A generator
    """
    while True:
        partial_file = file.read(size)
        if not partial_file:
            break
        yield partial_file


def write_out_parquet(df: pl.LazyFrame, outpath: Path, filename: str):
    """
    Writes out a dataframe to parquet format.

    Args:
        df (pl.DataFrame): A dataframe to write out
        outpath (Path): Path to the output folder
        filename (str): Name of the file
    """
    df.collect().write_parquet(outpath / f"{filename}.parquet")


def get_n_lines(file_path: Path) -> int:
    """
    Determines how many lines have to be processed
    depending on options and number of available lines.
    Checks that the number of lines is a multiple of 4.

    Args:
        file_path (string): Path to a fastq.gz file

    Returns:
        n_lines (int): Number of lines in the file
    """
    print(f"Counting number of reads in file {file_path}")
    with gzip.open(file_path, "rt", encoding="utf-8", errors="ignore") as infile:
        n_lines = sum(bl.count("\n") for bl in blocks(infile))
    if n_lines % 4 != 0:
        sys.exit(
            f"{file_path}'s number of lines is not a multiple of 4. The file "
            "might be corrupted.\n Exiting"
        )
    return n_lines


def get_read_paths(read1_path: str, read2_path: str) -> Tuple[list[Path], list[Path]]:
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
    _read1_paths = [Path(path) for path in read1_path.split(",")]
    _read2_paths = [Path(path) for path in read2_path.split(",")]
    if len(_read1_paths) != len(_read2_paths):
        sys.exit(
            f"Unequal number of read1 ({len(_read1_paths)}) and read2({len(_read2_paths)}) files provided"
            "\n Exiting"
        )
    all_files = _read1_paths + _read2_paths
    for file_path in all_files:
        if os.path.isfile(file_path):
            if os.access(file_path, os.R_OK):
                continue
            else:
                sys.exit(f"{file_path} is not readable. Exiting")
        else:
            sys.exit(f"{file_path} does not exist. Exiting")

    return (_read1_paths, _read2_paths)


def get_csv_reader_from_path(filename: str, sep: str = "\t"):
    """
    Returns a csv_reader object for a file weather it's a flat file or compressed.

    Args:
        filename: str

    Returns:
        csv_reader: The csv_reader for the file
    """
    if filename.endswith(".gz"):
        file_handle = gzip.open(filename, mode="rt")
        csv_reader = csv.reader(file_handle, delimiter=sep)
    else:
        file_handle = open(filename, encoding="UTF-8")
        csv_reader = csv.reader(file_handle, delimiter=sep)
    return csv_reader


def check_file(file_str: str) -> Path:
    """Check that a file exists and is readable.

    Args:
        file_str (str): Path to the file as a string

    Returns:
        Path: Path to the file
    """
    file_path = Path(file_str)
    if not file_path.exists:
        raise FileNotFoundError(f"This file {file_path} does not exist")

    if not access(file_path, R_OK):
        raise FileNotFoundError(f"This file {file_path} is not accessible")
    return file_path


def write_unmapped(unmapped_df: pl.LazyFrame, outfolder: Path, filename: str):
    """
    Writes a list of top unmapped sequences

    Args:
        merged_no_match (Counter): Counter of unmapped sequences
        top_unknowns (int): Number of unmapped sequences to output
        outfolder (string): Path of the output folder
        filename (string): Name of the output file
    """

    unmapped_df.collect().write_csv(file=outfolder / f"{filename}.csv")


def load_report_template() -> dict:
    """Load json template for the report

    Returns:
        dict: Dict for the report
    """
    with open(file=JSON_REPORT_PATH, encoding="utf-8") as json_in:
        try:
            report_dict = json.load(json_in)
        except json.JSONDecodeError:
            sys.exit(
                f"Json report template at {JSON_REPORT_PATH} is not a valid json file."
            )
    return report_dict


def create_report(
    total_reads: int,
    unmapped: pl.DataFrame,
    version: str,
    start_time,
    umis_corrected: int,
    bcs_corrected: int,
    bad_cells,
    r1_too_short: int,
    r2_too_short: int,
    args: Namespace,
    chemistry_def,
    maximum_distance: int,
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
    total_unmapped = unmapped[COUNT_COLUMN][0]
    total_too_short = r1_too_short + r2_too_short
    total_mapped = total_reads - total_unmapped - total_too_short

    too_short_perc = round((total_too_short / total_reads) * 100)
    mapped_perc = round((total_mapped / total_reads) * 100)
    unmapped_perc = round((total_unmapped / total_reads) * 100)

    report_data = load_report_template()
    report_data["Date"] = datetime.datetime.today().strftime("%Y-%m-%d")
    report_data["Running time"] = secondsToText.secondsToText(time.time() - start_time)
    report_data["CITE-seq-Count Version"] = version
    report_data["Reads processed"] = int(total_reads)
    report_data["Percentage mapped"] = mapped_perc
    report_data["Percentage unmapped"] = unmapped_perc
    report_data["Percentage too short"] = too_short_perc
    report_data["r1_too_short"] = r1_too_short
    report_data["r2_too_short"] = r2_too_short
    report_data["Uncorrected cells"] = len(bad_cells)
    report_data["Correction"]["Cell barcodes collapsing threshold"] = args.bc_threshold
    report_data["Correction"]["Cell barcodes corrected"] = bcs_corrected
    report_data["Correction"]["UMI collapsing threshold"] = args.umi_threshold
    report_data["Correction"]["UMIs corrected"] = umis_corrected
    report_data["Run parameters"]["Read1_paths"] = args.read1_path
    report_data["Run parameters"]["Read2_paths"] = args.read2_path
    report_data["Run parameters"]["Cell barcode"][
        "First position"
    ] = chemistry_def.cell_barcode_start
    report_data["Run parameters"]["Cell barcode"][
        "Last position"
    ] = chemistry_def.cell_barcode_end
    report_data["Run parameters"]["UMI barcode"][
        "First position"
    ] = chemistry_def.umi_barcode_start
    report_data["Run parameters"]["UMI barcode"][
        "Last position"
    ] = chemistry_def.umi_barcode_end
    report_data["Expected cells"] = args.expected_barcodes
    report_data["Tags max errors"] = maximum_distance
    report_data["Start trim"] = chemistry_def.r2_trim_start

    with open(
        os.path.join(args.outfolder, "run_report.yaml"), "w", encoding="utf-8"
    ) as report_file:
        yaml.dump(report_data, report_file, default_flow_style=False, sort_keys=False)


def create_mtx_df(
    main_df: pl.LazyFrame, tags_df: pl.DataFrame, subset_df: pl.LazyFrame
) -> tuple[pl.DataFrame, pl.DataFrame, pl.DataFrame]:
    """Create the MTX dataframe, indexed barcodes and indexed features from the different dataframes

    Args:
        main_df (pl.DataFrame): Bridge between barcodes and features. Holds UMIs
        tags_df (pl.DataFrame): Features and their sequences
        subset_df (pl.DataFrame): Subset of barcodes to use

    Returns:
        tuple[pl.DataFrame, pl.DataFrame, pl.DataFrame]: MTX dataframe, indexed barcodes, indexed features
    """
    tags_indexed = (
        pl.concat(
            [
                tags_df.lazy(),
                pl.LazyFrame(
                    {FEATURE_NAME_COLUMN: UNMAPPED_NAME, SEQUENCE_COLUMN: "UNKNOWN"}
                ),
            ]
        )
        .sort(pl.col(FEATURE_NAME_COLUMN))
        .with_row_count(offset=1, name=FEATURE_ID_COLUMN)
    )
    barcodes_indexed = subset_df.sort(pl.col(SUBSET_COLUMN)).with_row_count(
        offset=1, name=BARCODE_ID_COLUMN
    )
    mtx_df = (
        tags_indexed.join(main_df, on=FEATURE_NAME_COLUMN)
        .join(barcodes_indexed, left_on=BARCODE_COLUMN, right_on=SUBSET_COLUMN)
        .select([FEATURE_ID_COLUMN, BARCODE_ID_COLUMN, COUNT_COLUMN])
    )
    return mtx_df.collect(), tags_indexed.collect(), barcodes_indexed.collect()


def write_data_to_mtx(
    main_df: pl.LazyFrame,
    tags_df: pl.DataFrame,
    subset_df: pl.LazyFrame,
    data_type: str,
    outpath: str,
) -> None:
    """Write out the data to disk in MTX format

    Args:
        main_df (pl.DataFrame): Main df
        tags_df (pl.DataFrame): Features df
        subset_df (pl.DataFrame): Subsetted barcodes df
        data_type (str): umi or read
        outpath (str): Path to the output folder
    """
    mtx_df, tags_indexed, barcodes_indexed = create_mtx_df(main_df, tags_df, subset_df)
    data_path = Path(outpath) / f"{data_type}_count"
    data_path.mkdir(parents=True, exist_ok=True)
    mtx_df.write_csv(include_header=False, file=data_path / TEMP_MTX, separator=" ")
    # Write out the full MTX matrix
    first_line_mtx = (
        f"{tags_indexed.shape[0]} {barcodes_indexed.shape[0]} {mtx_df.shape[0]}\n"
    )
    with open(data_path / TEMP_MTX, "r") as mtx_in:
        mtx_main = mtx_in.read()
        final_mtx = MTX_HEADER + first_line_mtx + mtx_main
    with open(data_path / TEMP_MTX, "w") as mtx_out:
        mtx_out.write(final_mtx)
    with open(data_path / TEMP_MTX, "rb") as mtx_in:
        with gzip.open(data_path / MATRIX_MTX, "wb") as mtx_gz:
            shutil.copyfileobj(mtx_in, mtx_gz)
    os.remove(data_path / TEMP_MTX)
    # Write ouf features and barcodes
    tags_indexed.sort(FEATURE_ID_COLUMN).select(FEATURE_NAME_COLUMN).write_csv(
        file=data_path / "features.csv", include_header=False
    )
    compress_file(data_path / "features.csv", data_path / FEATURES_MTX)
    barcodes_indexed.sort(BARCODE_ID_COLUMN).select(SUBSET_COLUMN).write_csv(
        file=data_path / "barcodes.csv", include_header=False
    )
    compress_file(data_path / "barcodes.csv", data_path / BARCODE_MTX)


def write_mapping_input(
    args: Namespace,
    read1_paths: list[Path],
    read2_paths: list[Path],
    r2_min_length: int,
    chemistry_def,
):
    """
    Writes all reads to one CSV to be used.

    Args:
        args(argparse): All parsed arguments.
        read1_paths (list): List of R1 fastq.gz paths.
        read2_paths (list): List of R2 fastq.gz paths.
        r2_min_length (int):  Minimum length of read2 sequences.
        chemistry_def (namedtuple): Hols all the information about the chemistry definition.
        parsed_tags (list): List of namedtuple tags.
        maximum_distance (int): Maximum hamming distance for mapping.
    """
    print("Writing reads to disk")

    temp_path = os.path.abspath(args.temp_path)
    r1_too_short = 0
    r2_too_short = 0
    total_reads = 0
    total_reads_written = 0

    barcode_slice = slice(
        chemistry_def.cell_barcode_start - 1, chemistry_def.cell_barcode_end
    )
    umi_slice = slice(
        chemistry_def.umi_barcode_start - 1, chemistry_def.umi_barcode_end
    )

    temp_file = tempfile.NamedTemporaryFile(
        "w", dir=temp_path, suffix="_csc", delete=False
    )
    temp_file_path = Path(temp_file.name)
    reads_written = 0

    for read1_path, read2_path in zip(read1_paths, read2_paths):
        print(f"Reading reads from files: {read1_path}, {read2_path}")
        with gzip.open(read1_path, "rt") as textfile1, gzip.open(
            read2_path, "rt"
        ) as textfile2:
            secondlines = islice(zip(textfile1, textfile2), 1, None, 4)

            for read1, read2 in secondlines:
                total_reads += 1

                read1 = read1.strip()
                if len(read1) < chemistry_def.umi_barcode_end:
                    r1_too_short += 1
                    # The entire read is skipped
                    continue
                if len(read2) < r2_min_length:
                    r2_too_short += 1
                    # The entire read is skipped
                    continue

                read1_sliced = read1[
                    chemistry_def.cell_barcode_start - 1 : chemistry_def.umi_barcode_end
                ]

                read2_sliced = read2[
                    chemistry_def.r2_trim_start : (
                        r2_min_length + chemistry_def.r2_trim_start
                    )
                ]
                temp_file.write(
                    "{},{},{}\n".format(
                        read1_sliced[barcode_slice],
                        read1_sliced[umi_slice],
                        read2_sliced,
                    )
                )

                reads_written += 1
                total_reads_written += 1

    temp_file.close()

    return (
        temp_file_path,
        r1_too_short,
        r2_too_short,
        total_reads,
    )
