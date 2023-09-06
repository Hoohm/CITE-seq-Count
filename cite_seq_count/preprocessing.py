"""Sets of functions to preprocess the data"""


import gzip
import sys
from itertools import combinations, islice
import Levenshtein


from cite_seq_count.io import get_n_lines, check_file
import polars as pl
from pandas import read_csv

REQUIRED_TAGS_HEADER = ["sequence", "feature_name"]
REQUIRED_CELLS_REF_HEADER = ["reference"]
OPTIONAL_CELLS_REF_HEADER = ["translation"]
FEATURE_NAME = "feature_name"
SEQUENCE = "sequence"
REQUIRED_TAGS_HEADER = [FEATURE_NAME, SEQUENCE]
STRIP_CHARS = '"0123456789- \t\n'


def parse_cell_list_csv(filename: str, barcode_length: int) -> pl.DataFrame:
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
    cells_pl = pl.read_csv(file_path)
    barcode_pattern = rf"^[ATGC]{{{barcode_length}}}"

    header = cells_pl.columns
    set_dif = set(REQUIRED_CELLS_REF_HEADER) - set(header)
    if len(set_dif) != 0:
        set_diff_string = ",".join(list(set_dif))
        raise SystemExit(f"The header is missing {set_diff_string}. Exiting")
    if OPTIONAL_CELLS_REF_HEADER in header:
        with_translation = True
    else:
        with_translation = False
    # Prepare and validate cells_pl
    if with_translation:
        cells_pl = cells_pl.with_columns(
            reference=pl.col("reference").str.strip(STRIP_CHARS),
            translation=pl.col("translation").str.strip(STRIP_CHARS),
        )

    else:
        cells_pl = cells_pl.with_columns(
            reference=pl.col("reference").str.strip(STRIP_CHARS),
        )

    check_sequence_pattern(
        df=cells_pl,
        pattern=barcode_pattern,
        column_name="reference",
        file_type="Cell reference",
        expected_pattern="ATGC",
    )

    if with_translation:
        check_sequence_pattern(
            df=cells_pl,
            pattern=barcode_pattern,
            column_name="translation",
            file_type="Cell reference",
            expected_pattern="ATGC",
        )

    return cells_pl


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
        set_diff_str = " AND ".join(set_diff)
        raise SystemExit(
            f"The header of the tags file at {file_path} is missing the following header(s) {set_diff_str}"
        )
    for column in REQUIRED_TAGS_HEADER:
        if not data_pl.filter(pl.col([column]) is None).is_empty():
            raise SystemExit(
                f"Column {column} is missing a value. Please fix the CSV file."
            )
    check_sequence_pattern(
        df=data_pl,
        pattern=atgc_test,
        column_name=SEQUENCE,
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
            regex_test.filter(pl.col("regex") == False)
            .get_column(column_name)
            .to_list()
        )
        sequences_str = "\n".join(sequences)
        raise SystemExit(
            f"Some sequences in the {file_type} file is not only composed of the proper pattern {expected_pattern}.\nHere are the sequences{sequences_str}"
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

    # Check if the distance is big enoughbetween tags
    offending_pairs = []
    for tag_a, tag_b in combinations(tags_pl["sequence"], 2):
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
    longest_tag_len = max(tags_pl["sequence"].str.n_chars())

    return longest_tag_len


def get_read_length(filename):
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
                    "[ERROR] Sequence length in {} is not consistent. Please, trim all "
                    "sequences at the same length.\n"
                    "Exiting the application.\n".format(filename)
                )
    return read_length


def translate_barcodes(cell_set, translation_dict):
    """Translate a list of barcode using a mapping translation
    Args:
        cell_set (set): A set of barcodes
        translation_dict (dict): A dict providing a simple key value translation

    Returns:
        translated_barcodes (set): A set of translated barcodes
    """

    translated_barcodes = set()
    for translated_barcode in translation_dict.keys():
        if translation_dict[translated_barcode] in cell_set:
            translated_barcodes.add(translated_barcode)
    return translated_barcodes


def check_barcodes_lengths(read1_length, cb_first, cb_last, umi_first, umi_last):
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
            "[WARNING] Read1 length is {}bp but you are using {}bp for Cell "
            "and UMI barcodes combined.\nThis might lead to wrong cell "
            "attribution and skewed umi counts.\n".format(
                read1_length, barcode_umi_length
            )
        )


def pre_run_checks(read1_paths, chemistry_def, longest_tag_len, args):
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


def get_filtered_list(args, chemistry, translation_dict):
    """
    Determines what mode to use for cell barcode correction.
    Args:
        args(argparse): All arguments

    Returns:
        set if we have a filtered list
        None if we want correction and we have not a list
        False if we deactivation filtering
    """
    if args.filtered_cells:
        filtered_set = parse_filtered_list_csv(
            args.filtered_cells,
            (chemistry.cell_barcode_end - chemistry.cell_barcode_start),
        )
        # Do we need to translate the list?
        if args.translation_list:
            # get the translation
            translated_set = translate_barcodes(
                cell_set=filtered_set, translation_dict=translation_dict
            )
            return translated_set
        return filtered_set
    else:
        return None
