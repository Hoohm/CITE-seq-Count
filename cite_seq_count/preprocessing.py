"""Sets of functions to preprocess the data"""


import gzip
import sys
from itertools import combinations, islice
from collections import Counter
from argparse import ArgumentParser
from pathlib import Path
import Levenshtein
import polars as pl
import umi_tools.whitelist_methods as whitelist_method
from cite_seq_count.io import get_n_lines, check_file
from cite_seq_count.chemistry import Chemistry


# REQUIRED_TAGS_HEADER = ["sequence", "feature_name"]
REQUIRED_CELLS_REF_HEADER = ["reference"]
OPTIONAL_CELLS_REF_HEADER = ["translation"]
# Polars column names
# Tags input
FEATURE_NAME_COLUMN = "feature_name"
SEQUENCE_COLUMN = "sequence"
REQUIRED_TAGS_HEADER = [FEATURE_NAME_COLUMN, SEQUENCE_COLUMN]
# Reads input
BARCODE_COLUMN = "barcode"
CORRECTED_BARCODE_COLUMN = "corrected_barcode"
UMI_COLUMN = "umi"
R2_COLUMN = "r2"
# Barcode input
REFERENCE_COLUMN = "reference"
TRANSLATION_COLUMN = "translation"
WHITELIST_COLUMN = "whitelist"
STRIP_CHARS = '"0123456789- \t\n'

UNMAPPED_NAME = "unmapped"


def parse_barcode_reference(
    filename: str, barcode_length: int, required_header: str
) -> pl.DataFrame:
    """Reads white-listed barcodes from a CSV file.

    The function accepts plain barcodes or even 10X style barcodes with the
    `-1` at the end of each barcode.

    Args:
        filename (str): translation_list barcode file.
        barcode_length (int): Length of the expected barcodes.

    Returns:
        set: The set of white-listed barcodes.

    """
    file_path = check_file(filename)
    barcodes_pl = pl.read_csv(file_path.absolute())
    barcode_pattern = rf"^[ATGC]{{{barcode_length}}}"

    header = barcodes_pl.columns
    set_dif = set(required_header) - set(header)
    if len(set_dif) != 0:
        set_diff_string = ",".join(list(set_dif))
        raise SystemExit(f"The header is missing {set_diff_string}. Exiting")
    if OPTIONAL_CELLS_REF_HEADER in header:
        with_translation = True
    else:
        with_translation = False
    # Prepare and validate barcodes_pl
    if with_translation:
        barcodes_pl = barcodes_pl.with_columns(
            reference=pl.col(REFERENCE_COLUMN).str.strip_chars(STRIP_CHARS),
            translation=pl.col(TRANSLATION_COLUMN).str.strip_chars(STRIP_CHARS),
        )

    else:
        barcodes_pl = barcodes_pl.with_columns(
            reference=pl.col(REFERENCE_COLUMN).str.strip_chars(STRIP_CHARS),
        )

    check_sequence_pattern(
        df=barcodes_pl,
        pattern=barcode_pattern,
        column_name=REFERENCE_COLUMN,
        file_type="Barcode reference",
        expected_pattern="ATGC",
    )

    if with_translation:
        check_sequence_pattern(
            df=barcodes_pl,
            pattern=barcode_pattern,
            column_name=TRANSLATION_COLUMN,
            file_type="Barcode reference",
            expected_pattern="ATGC",
        )

    return barcodes_pl


def parse_tags_csv(file_name: str) -> pl.DataFrame:
    """Reads the TAGs from a CSV file. Checks that the header contains
    necessary strings and if sequences are made of ATGC

    The expected file format has a header with "sequence" and "feature_name".
    Order doesn't matter.
    e.g. file content
        sequence,feature_name
        GTCAACTCTTTAGCG,Hashtag_1
        TGATGGCCTATTGGG,Hashtag_2
        TTCCGCCTCTCTTTG,Hashtag_3

    Args:
        file_name (str): file path as a tring

    Returns:
        pl.DataFrame: polars dataframe with the csv content
    """

    atgc_test = "^[ATGC]{1,}$"

    file_path = check_file(file_str=file_name)

    data_pl = pl.read_csv(file_path)
    file_header = set(data_pl.columns)
    set_diff = set(REQUIRED_TAGS_HEADER).difference(file_header)
    if len(set_diff) != 0:
        print(set_diff)
        print(len(set_diff))
        set_diff_str = " AND ".join(set_diff)
        raise SystemExit(
            f"The header of the tags file at {file_path}"
            f"is missing the following header(s) {set_diff_str}"
        )
    for column in REQUIRED_TAGS_HEADER:
        if not data_pl.filter(pl.col([column]) is None).is_empty():
            raise SystemExit(
                f"Column {column} is missing a value. Please fix the CSV file."
            )
    check_sequence_pattern(
        df=data_pl,
        pattern=atgc_test,
        column_name=SEQUENCE_COLUMN,
        file_type="tags",
        expected_pattern="ATGC",
    )

    return data_pl


def check_sequence_pattern(
    df: pl.DataFrame,
    pattern: str,
    column_name: str,
    file_type: str,
    expected_pattern: str,
) -> None:
    """Check that a column of a polars df matches a given pattern and exit if not

    Args:
        df (pl.DataFrame): Df holding the info to be tested
        pattern (str): Regex pattern to be tested
        column_name (str): Which column to test
        file_type (str): File type for the error raised
        expected_pattern (str): Human readable pattern to be raised

    Raises:
        SystemExit: Exists if some patterns don't match
    """
    regex_test = df.with_columns(
        pl.col(column_name).str.contains(pattern).alias("regex")
    )
    if not regex_test.select(pl.col("regex").all()).get_column("regex").item():
        sequences = (
            regex_test.filter(pl.col("regex") is False)
            .get_column(column_name)
            .to_list()
        )
        sequences_str = "\n".join(sequences)
        raise SystemExit(
            f"Some sequences in the {file_type} file is not only composed"
            f"of the proper pattern {expected_pattern}.\n"
            f"Here are the sequences{sequences_str}"
        )


def check_tags(tags_pl: pl.DataFrame, maximum_distance: int) -> int:
    """Evaluates the distance between the TAGs based on the `maximum distance`
    argument provided.

    The output will have the keys sorted by TAG length (longer first). This
    way, longer barcodes will be evaluated first.
    Adds unmapped category as well.

    Args:
        tags (dict): A dictionary with TAG sequences as keys and names as values.
        maximum_distance (int): The minimum Levenshtein distance allowed
            between two TAGs.

    Returns:
        int: the length of the longest TAG

    """
    # TODO: Decide to keep or delete.
    # # Check that all tags are the same length
    # if tags_pl[SEQUENCE_COLUMN].str.lengths().unique().shape[0] != 1:
    #     SystemExit(
    #         "Tag sequences have different lengths. Version 2 can only run with one
    #          length. Please use an older version"
    #     )
    # Check if the distance is big enough between tags
    offending_pairs = []
    for tag_a, tag_b in combinations(tags_pl[SEQUENCE_COLUMN], 2):
        # pylint: disable=no-member
        distance = Levenshtein.distance(tag_a, tag_b)
        if distance <= (maximum_distance - 1):
            offending_pairs.append([tag_a, tag_b, distance])
    # If offending pairs are found, print them all.
    if offending_pairs:
        print(
            "[ERROR] Minimum Levenshtein distance of TAGs barcode is less "
            "than given threshold.\n"
            "Please use a smaller distance.\n\n"
            "Offending case(s):\n"
        )
        for pair in offending_pairs:
            print(f"\t{pair[0]}\n\t{pair[1]}\n\tDistance = {pair[2]}\n")
        sys.exit("Exiting the application.\n")
    longest_tag_len = max(tags_pl[SEQUENCE_COLUMN].str.n_chars())

    return longest_tag_len


def get_read_length(filename: Path):
    """Check wether SEQUENCE lengths are consistent in
    the first 1000 reads from a FASTQ file and return
    the length.

    Args:
        filename (str): FASTQ file.

    Returns:
        int: The file's SEQUENCE length.

    """
    with gzip.open(filename, "r") as fastq_file:
        secondlines = islice(fastq_file, 1, 1000, 4)
        temp_length = len(next(secondlines).rstrip())
        for sequence in secondlines:
            read_length = len(sequence.rstrip())
            if temp_length != read_length:
                sys.exit(
                    f"[ERROR] Sequence length in {filename} is not consistent."
                    f" Please, trim all sequences at the same length.\n"
                    f"Exiting the application.\n"
                )
    return read_length


def check_barcodes_lengths(
    read1_length: int, cb_first: int, cb_last: int, umi_first: int, umi_last: int
):
    """Check Read1 length against CELL and UMI barcodes length.

    Args:
        read1_length (int): Read1 length.
        cb_first (int): Barcode first base position for Read1.
        cb_last (int): Barcode last base position for Read1.
        umi_first (int): UMI first base position for Read1.
        umi_last (int): UMI last base position for Read1.
    """
    barcode_length = cb_last - cb_first + 1
    umi_length = umi_last - umi_first + 1
    barcode_umi_length = barcode_length + umi_length

    if barcode_umi_length > read1_length:
        raise SystemExit(
            "[ERROR] Read1 length is shorter than the option you are using for "
            "Cell and UMI barcodes length. Please, check your options and rerun.\n\n"
            "Exiting the application.\n"
        )
    elif barcode_umi_length < read1_length:
        print(
            f"[WARNING] Read1 length is {read1_length}bp"
            f"but you are using {barcode_umi_length}bp for Cell "
            f"and UMI barcodes combined.\nThis might lead to wrong cell "
            f"attribution and skewed umi counts.\n"
        )


def pre_run_checks(
    read1_paths: list[Path],
    chemistry_def: dict,
    longest_tag_len: int,
    args: ArgumentParser,
):
    """Checks that the chemistry is properly set and defines how many reads to process

    Args:
        read1_paths (list): List of paths
        chemistry_def (Chemistry): Chemistry definition
        longest_tag_len (int): Longest tag sequence
        args (argparse): List of arguments

    Returns:
        n_reads (int): Number of reads to run on
        R2_min_length (int): Min R2 length to check if reads are too short
        maximum_distance (int): Maximum error rate allowed for mapping tags

    """
    read1_lengths = []
    read2_lengths = []
    total_reads = 0

    for read1_path in read1_paths:
        n_lines = get_n_lines(read1_path)
        total_reads += n_lines / 4
        # Get reads length. So far, there is no validation for Read2.
        read1_lengths.append(get_read_length(read1_path))

        # Check Read1 length against CELL and UMI barcodes length.
        check_barcodes_lengths(
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

    number_of_samples = len(read1_paths)

    # Print a statement if multiple files are run.
    if number_of_samples != 1:
        print(f"Detected {number_of_samples} pairs of files to run on.")

    if args.sliding_window:
        r2_min_length = read2_lengths[0]
        maximum_distance = 0
    else:
        r2_min_length = longest_tag_len
        maximum_distance = args.max_error
    return n_reads, r2_min_length, maximum_distance


def split_data_input(mapping_input_path: Path):
    input_df = (
        pl.read_csv(
            mapping_input_path,
            has_header=False,
            new_columns=[BARCODE_COLUMN, UMI_COLUMN, R2_COLUMN],
        )
        .group_by([BARCODE_COLUMN, UMI_COLUMN, R2_COLUMN])
        .agg(pl.count())
    )

    barcodes_df = (
        input_df.select([BARCODE_COLUMN, "count"])
        .group_by(BARCODE_COLUMN)
        .agg(pl.sum("count"))
    )
    r2_df = input_df.select(R2_COLUMN).unique()

    return input_df, barcodes_df, r2_df


def get_barcode_subset(
    barcode_whitelist: Path,
    expected_barcodes: int,
    chemistry: Chemistry,
    barcode_reference: pl.DataFrame | None,
    barcodes_df: pl.DataFrame,
):
    """
    Generate the barcode list used for barcode correction and subsetting
    """
    enable_barcode_correction = True
    if barcode_whitelist:
        barcode_subset = parse_barcode_reference(
            filename=expected_barcodes,
            barcode_length=(chemistry.cell_barcode_end - chemistry.cell_barcode_start),
            required_header=WHITELIST_COLUMN,
        )
    else:
        n_barcodes = barcode_whitelist
        if barcode_reference is not None:
            barcode_subset = (
                barcodes_df.filter(
                    pl.col(BARCODE_COLUMN).str.is_in(
                        barcode_reference[REFERENCE_COLUMN]
                    )
                )
                .group_by(BARCODE_COLUMN)
                .agg(pl.count())
                .sort("count", descending=True)
                .head(n_barcodes * 1.2)
                .drop("count")
                .rename({SEQUENCE_COLUMN: WHITELIST_COLUMN})
            )
        else:
            raw_barcodes_dict = (
                barcodes_df.filter(~pl.col(BARCODE_COLUMN).str.contains("N"))
                .group_by(BARCODE_COLUMN)
                .agg(pl.count())
                .sort("count", descending=True)
            ).to_dict()
            barcode_counter = Counter(
                zip(raw_barcodes_dict[BARCODE_COLUMN], raw_barcodes_dict["count"])
            )
            true_barcodes = whitelist_method.getKneeEstimateDistance(
                cell_barcode_counts=barcode_counter, cell_number=n_barcodes
            )
            barcode_subset = pl.DataFrame(
                true_barcodes, schema={WHITELIST_COLUMN: pl.Utf8, "counts": pl.UInt32}
            ).drop("counts")

    if n_barcodes > barcode_subset.shape[0]:
        print(
            f"Number of expected cells, {n_barcodes}, is higher "
            f"than number of cells found {barcode_subset.shape[0]}.\nNot performing "
            f"cell barcode correction"
        )
        enable_barcode_correction = False
    return barcode_subset, enable_barcode_correction
