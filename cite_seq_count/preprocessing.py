"""Sets of functions to preprocess the data"""


import gzip
import sys
from itertools import combinations, islice
from collections import Counter
from argparse import Namespace
from pathlib import Path
import Levenshtein
import polars as pl
import umi_tools.whitelist_methods as whitelist_method
from cite_seq_count.io import get_n_lines, check_file
from cite_seq_count.constants import (
    SEQUENCE_COLUMN,
    R2_COLUMN,
    BARCODE_COLUMN,
    UMI_COLUMN,
    REQUIRED_TAGS_HEADER,
    REFERENCE_COLUMN,
    TRANSLATION_COLUMN,
    OPTIONAL_CELLS_REF_HEADER,
    STRIP_CHARS,
    SUBSET_COLUMN,
    COUNT_COLUMN,
)


def check_equi_length(df: pl.DataFrame, column_name: str):
    """Check that all the barcodes in the specified column of a polars DataFrame are the same length.

    Args:
        df (pl.DataFrame): The DataFrame containing the barcodes.
        column_name (str): The name of the column containing the barcodes.

    Raises:
        ValueError: If the barcodes have different lengths.

    """
    barcode_lengths = df[column_name].str.len_chars().unique()
    if len(barcode_lengths) > 1:
        raise ValueError(f"Barcodes in {column_name} column have different lengths.")


def parse_barcode_file(
    filename: str, barcode_length: int, required_header: list
) -> pl.DataFrame:
    """Reads reference barcodes from a CSV file.

    The function accepts plain barcodes or even 10X style barcodes with the
    `-1` at the end of each barcode.

    Args:
        filename (str): translation_list barcode file.
        barcode_length (int): Length of the expected barcodes.

    Returns:
        set: The set of reference barcodes.

    """
    file_path = check_file(filename)
    barcodes_df = pl.read_csv(file_path.absolute())
    barcode_pattern = "^[ATGC]{1,}$"

    header = barcodes_df.columns
    set_dif = set(required_header) - set(header)
    if len(set_dif) != 0:
        set_diff_string = ",".join(list(set_dif))
        raise SystemExit(f"The header is missing {set_diff_string}. Exiting")
    if OPTIONAL_CELLS_REF_HEADER in header:
        with_translation = True
    else:
        with_translation = False
    # Prepare and validate barcodes_df
    if with_translation:
        barcodes_df = barcodes_df.with_columns(
            reference=pl.col(REFERENCE_COLUMN).str.strip_chars(STRIP_CHARS),
            translation=pl.col(TRANSLATION_COLUMN).str.strip_chars(STRIP_CHARS),
        )

    else:
        barcodes_df = barcodes_df.with_columns(
            reference=pl.col(REFERENCE_COLUMN).str.strip_chars(STRIP_CHARS),
        )

    check_sequence_pattern(
        df=barcodes_df,
        pattern=barcode_pattern,
        column_name=REFERENCE_COLUMN,
        file_type="Barcode reference",
        expected_pattern="ATGC",
        filename=filename,
    )

    check_equi_length(df=barcodes_df, column_name=REFERENCE_COLUMN)

    if with_translation:
        check_sequence_pattern(
            df=barcodes_df,
            pattern=barcode_pattern,
            column_name=TRANSLATION_COLUMN,
            file_type="Barcode reference",
            expected_pattern="ATGC",
            filename=filename,
        )
        check_equi_length(df=barcodes_df, column_name=TRANSLATION_COLUMN)

    return barcodes_df


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
        filename=file_name,
    )
    return data_pl


def check_sequence_pattern(
    df: pl.DataFrame,
    pattern: str,
    column_name: str,
    file_type: str,
    expected_pattern: str,
    filename: str,
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
            f"Some sequences in the {file_type} file is not only composed "
            f"of {expected_pattern} in the column: {column_name}. "
            f"Here are the sequences: {sequences_str}. Filepath: {filename}"
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


def get_read_length(filename: Path) -> int:
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
        read_length = 0
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
    chemistry_def,
    longest_tag_len: int,
    arguments: Namespace,
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
    if arguments.first_n < float("inf"):
        n_reads = arguments.first_n
    else:
        n_reads = total_reads

    number_of_samples = len(read1_paths)

    # Print a statement if multiple files are run.
    if number_of_samples != 1:
        print(f"Detected {number_of_samples} pairs of files to run on.")

    # if arguments.sliding_window:
    #     r2_min_length = read2_lengths[0]
    #     maximum_distance = 0
    # else:
    r2_min_length = longest_tag_len
    maximum_distance = arguments.max_error
    return n_reads, r2_min_length, maximum_distance


def split_data_input(
    mapping_input_path: Path,
) -> tuple[pl.DataFrame, pl.DataFrame, pl.DataFrame]:
    """Read in all the input data and split it into three dataframes.

    Reduce the size of the data by grouping on barcodes, umis and sequences.
    The function splits the data into three main dataframes.
    1. main_df is the dataframe holding the links between the other dfs
    2. barcodes_df holds only the barcode information for barcode correction
    3. r2_df only holds the sequences for mapping

    Args:
        mapping_input_path (Path): Path to the csv file containing all the input data

    Returns:
        tuple[pl.DataFrame, pl.DataFrame, pl.DataFrame]: Three dfs described above
    """
    main_df = (
        pl.read_csv(
            mapping_input_path,
            has_header=False,
            new_columns=[BARCODE_COLUMN, UMI_COLUMN, R2_COLUMN],
        )
        .group_by([BARCODE_COLUMN, UMI_COLUMN, R2_COLUMN])
        .agg(pl.count())
    )

    barcodes_df = (
        main_df.select([BARCODE_COLUMN, COUNT_COLUMN])
        .group_by(BARCODE_COLUMN)
        .agg(pl.sum(COUNT_COLUMN))
    )
    r2_df = main_df.select(R2_COLUMN).unique()

    return main_df, barcodes_df, r2_df


def get_barcode_subset(
    barcode_reference: pl.DataFrame | None,
    n_barcodes: int,
    chemistry,
    barcode_subset: pl.DataFrame | None,
    barcodes_df: pl.DataFrame,
) -> tuple[pl.DataFrame, bool]:
    """Generate the barcode list used for barcode correction and subsetting

    Args:
        barcode_reference (Path): Barcode reference df (can hold translation column)
        expected_barcodes (int): Number of expected barcodes from user
        chemistry (Chemistry): Chemistry definition
        barcode_subset (pl.DataFrame | None): Specific subset given by the user
        barcodes_df (pl.DataFrame): Barcodes from the input data

    Returns:
        tuple[pl.DataFrame, bool]: Barcode subset, enable barcode correction
    """
    enable_barcode_correction = True
    # Subset: True
    if barcode_subset is None:
        # Subset: False, Reference: True
        if barcode_reference is not None:
            barcode_subset = (
                barcodes_df.filter(
                    pl.col(BARCODE_COLUMN).is_in(barcode_reference[REFERENCE_COLUMN])
                )
                .sort(COUNT_COLUMN, descending=True)
                .head(round(n_barcodes * 1.2))
                .drop(COUNT_COLUMN)
                .rename({BARCODE_COLUMN: SUBSET_COLUMN})
            )
        else:
            # Subset: False, Reference: False
            barcode_subset = find_knee_estimated_barcodes(barcodes_df=barcodes_df)

    if n_barcodes > barcode_subset.shape[0]:
        print(
            f"Number of expected cells, {n_barcodes}, is higher "
            f"than number of cells found {barcode_subset.shape[0]}.\nNot performing "
            f"cell barcode correction"
        )
        enable_barcode_correction = False
    return barcode_subset, enable_barcode_correction


def find_knee_estimated_barcodes(barcodes_df: pl.DataFrame) -> pl.DataFrame:
    """Find the subset of barcodes by the knee method

    Args:
        barcodes_df (pl.DataFrame): barcodes to use

    Returns:
        pl.DataFrame: Final list of barcodes
    """
    raw_barcodes_dict = barcodes_df.filter(
        ~pl.col(BARCODE_COLUMN).str.contains("N")
    ).sort("count", descending=True)
    barcode_counter = Counter()
    barcode_counts = dict(raw_barcodes_dict.iter_rows())  # type: ignore
    barcode_counter.update(barcode_counts)
    true_barcodes = whitelist_method.getKneeEstimateDistance(
        cell_barcode_counts=barcode_counter
    )
    barcode_subset = pl.DataFrame(true_barcodes, schema={SUBSET_COLUMN: pl.Utf8})
    return barcode_subset
