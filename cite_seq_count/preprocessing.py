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


def parse_whitelist_csv(csv_reader, barcode_length, collapsing_threshold):
    """Reads white-listed barcodes from a CSV file.

    The function accepts plain barcodes or even 10X style barcodes with the
    `-1` at the end of each barcode.

    Args:
        filename (str): Whitelist barcode file.
        barcode_length (int): Length of the expected barcodes.
        collapsing_threshold (int): Maximum distance to collapse cell barcodes.

    Returns:
        set: The set of white-listed barcodes.
        int: Collasping threshold

    """
    STRIP_CHARS = '"0123456789- \t\n'
    cell_pattern = regex.compile(r"[ATGC]{{{}}}".format(barcode_length))

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
    # collapsing_threshold=test_cell_distances(whitelist, collapsing_threshold)
    if len(whitelist) == 0:
        sys.exit(
            "Please check cell barcode indexes -cbs, -cbl because none of the given whitelist is valid."
        )
    return (set(whitelist), collapsing_threshold)


def test_cell_distances(whitelist, collapsing_threshold):
    """Tests cell barcode distances to validate provided cell barcode collapsing threshold
    
    Function needs the given whitelist as well as the threshold.
    If the value is too high, it will rerun until an acceptable value is found.
    
    Args:
        whitelist (set): Whitelist barcode set
        collapsing_threshold (int): Value of threshold

    Returns:
        collapsing_threshold (int): Valid threshold
    """
    ok = False
    while not ok:
        print(
            "Testing cell barcode collapsing threshold of {}".format(
                collapsing_threshold
            )
        )
        all_comb = combinations(whitelist, 2)
        for comb in all_comb:
            # pylint: disable=no-member
            if Levenshtein.hamming(comb[0], comb[1]) <= collapsing_threshold:
                collapsing_threshold -= 1
                print("Value is too high, reducing it by 1")
                break
        else:
            ok = True
    print("Using {} for cell barcode collapsing threshold".format(collapsing_threshold))
    return collapsing_threshold


def parse_tags_csv(filename):
    """Reads the TAGs from a CSV file.

    The expected file format (no header) is: TAG,TAG_NAME.
    e.g. file content
        GTCAACTCTTTAGCG,Hashtag_1
        TGATGGCCTATTGGG,Hashtag_2
        TTCCGCCTCTCTTTG,Hashtag_3

    Args:
        filename (str): TAGs file.

    Returns:
        dict: A dictionary containing the TAGs and their names.

    """
    with open(filename, mode="r") as csv_file:
        csv_reader = csv.reader(csv_file)
        tags = {}
        for row in csv_reader:
            tags[row[0].strip()] = row[1].strip()
    return tags


def check_tags(tags, maximum_distance):
    """Evaluates the distance between the TAGs based on the `maximum distance`
    argument provided.

    Additionally, it adds the barcode to the name of the TAG circumventing the
    need of having to share the mapping of the antibody and the barcode.
    
    The output will have the keys sorted by TAG length (longer first). This
    way, longer barcodes will be evaluated first.

    Args:
        tags (dict): A dictionary with the TAGs + TAG Names.
        maximum_distance (int): The maximum Levenshtein distance allowed
            between two TAGs.

    Returns:
        OrderedDict: An ordered dictionary containing the TAGs and
            their names in descendent order based on the length of the TAGs.
        int: the length of the longest TAG

    """
    ordered_tags = OrderedDict()
    longest_tag_len = 0
    for i, tag_seq in enumerate(sorted(tags, key=len, reverse=True)):
        ordered_tags[tags[tag_seq]] = {}
        ordered_tags[tags[tag_seq]]["id"] = i
        ordered_tags[tags[tag_seq]]["sequence"] = tag_seq
        if len(tag_seq) > longest_tag_len:
            longest_tag_len = len(tag_seq)

    ordered_tags["unmapped"] = {}
    ordered_tags["unmapped"]["id"] = i + 1
    ordered_tags["unmapped"]["sequence"] = "UNKNOWN"
    # If only one TAG is provided, then no distances to compare.
    if len(tags) == 1:
        ordered_tags["unmapped"] = {}
        ordered_tags["unmapped"]["id"] = 2
        return (ordered_tags, longest_tag_len)

    offending_pairs = []
    for a, b in combinations(tags.keys(), 2):
        # pylint: disable=no-member
        distance = Levenshtein.distance(a, b)
        if distance <= (maximum_distance - 1):
            offending_pairs.append([a, b, distance])
    DNA_pattern = regex.compile("^[ATGC]*$")
    for tag in tags:
        if not DNA_pattern.match(tag):
            print(
                "This tag {} is not only composed of ATGC bases.\nPlease check your tags file".format(
                    tag
                )
            )
            sys.exit("Exiting the application.\n")
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

    return (ordered_tags, longest_tag_len)


def sanitize_name(string):
    return string.replace("-", "_")


def convert_to_named_tuple(ordered_tags):
    # all_tags = namedtuple('all_tags', [sanitize_name(tag) for tag in ordered_tags.keys()])
    tag = namedtuple("tag", ["safe_name", "name", "sequence", "id"])
    tag_list = []
    for index, tag_name in enumerate(ordered_tags):
        tag_list.append(
            tag(
                safe_name=sanitize_name(tag_name),
                name=tag_name,
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


def get_chunk_strategy(read1_paths, read2_paths, chunk_size):
    pass


def check_barcodes_lengths(read1_length, cb_first, cb_last, umi_first, umi_last):
    """Check Read1 length against CELL and UMI barcodes length.

    Args:
        read1_length (int): Read1 length.
        cb_first (int): Barcode first base position for Read1.
        cb_last (int): Barcode last base position for Read1.
        umi_first (int): UMI first base position for Read1.
        umi_last (int): UMI last base position for Read1.

    Returns:
        slice: A `slice` object to extract the Barcode from the sequence string.
        slice: A `slice` object to extract the UMI from the sequence string.
        int: The Barcode + UMI length.

    """
    barcode_length = cb_last - cb_first + 1
    umi_length = umi_last - umi_first + 1
    barcode_umi_length = barcode_length + umi_length
    barcode_slice = slice(cb_first - 1, cb_last)
    umi_slice = slice(umi_first - 1, umi_last)

    if barcode_umi_length > read1_length:
        sys.exit(
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

    return (barcode_slice, umi_slice, barcode_umi_length)


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
    print("Counting number of reads")
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
            "\n Exiting".format(len(read1_path), len(read2_path))
        )
    return (_read1_path, _read2_path)
