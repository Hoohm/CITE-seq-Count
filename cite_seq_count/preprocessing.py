import csv
import gzip
import sys
import regex
import Levenshtein
import requests

from math import floor
from collections import OrderedDict
from collections import namedtuple
from itertools import combinations
from itertools import islice


def parse_whitelist_csv(filename, barcode_length):
    """Reads white-listed barcodes from a CSV file.

    The function accepts plain barcodes or even 10X style barcodes with the
    `-1` at the end of each barcode.

    Args:
        filename (str): Whitelist barcode file.
        barcode_length (int): Length of the expected barcodes.

    Returns:
        set: The set of white-listed barcodes.

    """
    STRIP_CHARS = '"0123456789- \t\n'
    cell_pattern = regex.compile(r"[ATGC]{{{}}}".format(barcode_length))

    if filename.endswith(".gz"):
        f = gzip.open(filename, mode="rt")
        csv_reader = csv.reader(f)
    else:
        f = open(filename, encoding="UTF-8")
        csv_reader = csv.reader(f)
    whitelist = [
        row[0].strip(STRIP_CHARS)
        for row in csv_reader
        if (len(row[0].strip(STRIP_CHARS)) == barcode_length)
    ]

    for cell_barcode in whitelist:
        if not cell_pattern.match(cell_barcode):
            sys.exit(
                "This barcode {} is not only composed of ATGC bases.".format(
                    cell_barcode
                )
            )
    if len(whitelist) == 0:
        sys.exit(
            "Please check cell barcode indexes -cbs, -cbl because none of the given whitelist is valid."
        )
    return set(whitelist)


def parse_tags_csv(filename):
    """Reads the TAGs from a CSV file. Checks if sequences are made of ATGC

    The expected file format has a header with "sequence" and "feature_name".
    Order doesn't matter.
    e.g. file content
        sequence,feature_name
        GTCAACTCTTTAGCG,Hashtag_1
        TGATGGCCTATTGGG,Hashtag_2
        TTCCGCCTCTCTTTG,Hashtag_3

    Args:
        filename (str): TAGs file.

    Returns:
        dict: A dictionary containing using sequences as keys and names as values.

    """
    REQUIRED_HEADER = ["sequence", "feature_name"]
    atgc_test = regex.compile("^[ATGC]{1,}$")
    with open(filename, mode="r") as csv_file:
        csv_reader = csv.reader(csv_file)
        tags = {}
        header = next(csv_reader)
        set_dif = set(REQUIRED_HEADER) - set(header)
        if len(set_dif) != 0:
            raise SystemExit(
                "The header is missing {}. Exiting".format(",".join(list(set_dif)))
            )
        sequence_id = header.index("sequence")
        feature_id = header.index("feature_name")
        for i, row in enumerate(csv_reader):
            sequence = row[sequence_id].strip()

            if not regex.match(atgc_test, sequence):
                raise SystemExit(
                    "Sequence {} on line {} is not only composed of ATGC. Exiting".format(
                        sequence, i
                    )
                )
            tags[sequence] = row[feature_id].strip()
    return tags


def check_tags(tags, maximum_distance):
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
        OrderedDict: An ordered dictionary containing the TAGs and
            their names in descendent order based on the length of the TAGs.
        int: the length of the longest TAG

    """
    tag = namedtuple("tag", ["name", "sequence", "id"])
    longest_tag_len = 0
    seq_list = []
    tag_list = []
    for i, tag_seq in enumerate(sorted(tags, key=len, reverse=True)):
        safe_name = sanitize_name(tags[tag_seq])

        # for index, tag_name in enumerate(ordered_tags):
        tag_list.append(tag(name=safe_name, sequence=tag_seq, id=i,))
        if len(tag_seq) > longest_tag_len:
            longest_tag_len = len(tag_seq)
        seq_list.append(tag_seq)
    tag_list.append(tag(name="unmapped", sequence="UNKNOWN", id=i + 1,))
    # If only one TAG is provided, then no distances to compare.
    if len(tags) == 1:
        return (tag_list, longest_tag_len)

    # Check if the distance is big enoughbetween tags
    offending_pairs = []
    for a, b in combinations(seq_list, 2):
        # pylint: disable=no-member
        distance = Levenshtein.distance(a, b)
        if distance <= (maximum_distance - 1):
            offending_pairs.append([a, b, distance])
    # If offending pairs are found, print them all.
    if offending_pairs:
        print(
            "[ERROR] Minimum Levenshtein distance of TAGs barcode is less "
            "than given threshold.\n"
            "Please use a smaller distance.\n\n"
            "Offending case(s):\n"
        )
        for pair in offending_pairs:
            print(
                "\t{tag1}\n\t{tag2}\n\tDistance = {distance}\n".format(
                    tag1=pair[0], tag2=pair[1], distance=pair[2]
                )
            )
        sys.exit("Exiting the application.\n")

    return (tag_list, longest_tag_len)


def sanitize_name(string):
    return string.replace("-", "_")


def convert_to_named_tuple(ordered_tags):
    # all_tags = namedtuple('all_tags', [sanitize_name(tag) for tag in ordered_tags.keys()])
    tag = namedtuple("tag", ["name", "sequence", "id"])
    tag_list = []
    for index, tag_name in enumerate(ordered_tags):
        tag_list.append(
            tag(
                name=ordered_tags[tag_name]["feature_name"],
                sequence=ordered_tags[tag_name]["sequence"],
                id=(index),
            )
        )
        # all_tags[index+1]=ordered_tags[tag_name]['sequence']
    return tag_list


def get_read_length(filename):
    """Check wether SEQUENCE lengths are consistent in a FASTQ file and return
    the length.

    Args:
        filename (str): FASTQ file.

    Returns:
        int: The file's SEQUENCE length.

    """
    with gzip.open(filename, "r") as fastq_file:
        secondlines = islice(fastq_file, 1, 1000, 4)
        # temp_length = len(next(secondlines).rstrip())
        for sequence in secondlines:
            read_length = len(sequence.rstrip())
            # if (temp_length != read_length):
            #     sys.exit(
            #         '[ERROR] Sequence length in {} is not consistent. Please, trim all '
            #         'sequences at the same length.\n'
            #         'Exiting the application.\n'.format(filename)
            #     )
    return read_length


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
    return (_read1_path, _read2_path)

