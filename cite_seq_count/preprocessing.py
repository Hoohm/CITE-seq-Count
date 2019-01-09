import csv
import gzip
import sys
import regex

import Levenshtein

from collections import OrderedDict
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
    with open(filename, mode='r') as csv_file:
        csv_reader = csv.reader(csv_file)
        whitelist = [row[0].strip(STRIP_CHARS) for row in csv_reader
                     if (len(row[0].strip(STRIP_CHARS)) == barcode_length)]
    return set(whitelist)


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
    with open(filename, mode='r') as csv_file:
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
        collections.OrderedDict: An ordered dictionary containing the TAGs and
            their names in descendent order based on the length of the TAGs.

    """
    ordered_tags = OrderedDict()
    for tag in sorted(tags, key=len, reverse=True):
        ordered_tags[tag] = tags[tag] + '-' + tag

    # If only one TAG is provided, then no distances to compare.
    if (len(tags) == 1):
        return(ordered_tags)
    
    offending_pairs = []
    for a, b in combinations(ordered_tags.keys(), 2):
        distance = Levenshtein.distance(a, b)
        if (distance <= maximum_distance):
            offending_pairs.append([a, b, distance])
    
    # If offending pairs are found, print them all.
    if offending_pairs:
        print(
            '[ERROR] Minimum Levenshtein distance of TAGs barcode is less '
            'than given threshold.\n'
            'Please use a smaller distance.\n\n'
            'Offending case(s):\n'
        )
        for pair in offending_pairs:
            print(
                '\t{tag1}\n\t{tag2}\n\tDistance = {distance}\n'
                .format(tag1=pair[0], tag2=pair[1], distance=pair[2])
            )
        sys.exit('Exiting the application.\n')
    
    return(ordered_tags)


def get_read_length(filename):
    """Check wether SEQUENCE lengths are consistent in a FASTQ file and return
    the length.

    Args:
        filename (str): FASTQ file.

    Returns:
        int: The file's SEQUENCE length.

    """
    with gzip.open(filename, 'r') as fastq_file:
        secondlines = islice(fastq_file, 1, 1000, 4)
        temp_length = len(next(secondlines).rstrip())
        for sequence in secondlines:
            read_length = len(sequence.rstrip())
            if (temp_length != read_length):
                sys.exit(
                    '[ERROR] Sequence length in {} is not consistent. Please, trim all '
                    'sequences at the same length.\n'
                    'Exiting the application.\n'.format(filename)
                )
    return(read_length)


def check_read_lengths(read1_length, cb_first, cb_last, umi_first, umi_last):
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
            '[ERROR] Read1 length is shorter than the option you are using for '
            'Cell and UMI barcodes length. Please, check your options and rerun.\n\n'
            'Exiting the application.\n'
        )
    elif barcode_umi_length < read1_length:
        print(
            '[WARNING] Read1 length is {}bp but you are using {}bp for Cell '
            'and UMI barcodes combined.\nThis might lead to wrong cell '
            'attribution and skewed umi counts.\n'
            .format(read1_length, barcode_umi_length)
        )
    
    return(barcode_slice, umi_slice, barcode_umi_length)

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
        if not b: break
        yield b


def get_n_lines(file_path, top_n):
    """
    Determines how many lines have to be processed
    depending on options and number of available lines.
    Checks that the number of lines is a multiple of 4.

    Args:
        file_path (string): Path to a fastq.gz file
        top_n (int): Number of reads to be used as defined by the user.

    Returns:
        n_lines (int): Number of lines to be used
    """
    print('Counting number of reads')
    with gzip.open(file_path, "rt",encoding="utf-8",errors='ignore') as f:
        n_lines = sum(bl.count("\n") for bl in blocks(f))
    if n_lines %4 !=0:
        sys.exit('{}\'s number of lines is not a multiple of 4. The file might be corrupted.\n Exiting')
    return(n_lines)


def generate_regex(tags, maximum_distance, error_type, legacy=False, max_poly_a=6, read2_length=98, user_regex=None):
    """Generate regex based ont he provided TAGs.

    Args:
        tags (dict): A dictionary with the TAGs + TAG Names.
        maximum_distance (int): The maximum Levenshtein distance allowed
            between two TAGs.
        legacy (bool): `True` if you use an earlier version of the kit that adds
            a T, C, or G at the end and you expect polyA tails in the data.
            Default is False.
        max_poly_a (int): Run length of A's expected for the polyA tail. Default
            is 6.
        read2_length (int): Length of Read2. Default is 98.
        user_regex (str): A regular expression to use for TAG matching. Default
            is None.

    Returns:
        regex.Pattern: An object that matches against any of the provided TAGs
            within the maximum distance provided.

    """
    # Get a list of the available TAGs.
    tag_keys = tags.keys()
    tags_set = set()
    for tag in tags:
        tags_set.add(len(tag))
    if len(tags_set) > 1:
        sys.exit(
            'Runnning CITE-seq-Count with tags of different lengths is not supported.\n'\
            'Please split your tags in different tag files and run them separately.')

    # Get the length of the longest TAG.
    longest_ab_tag = len(next(iter(tags)))

    if user_regex:
        # If more than one TAG is provided and their length is different, issue a
        # warning.
        if len(tag_keys) > 1:
            for i in range(1, len(tag_keys)):
                if len(tag_keys[i]) != len(tag_keys[i - 1]):
                    print(
                        '[WARNING] Different length TAGs have been provided while '
                        'you specified a custom Regex. An OR method is recommended '
                        'for this scenarios. No additional validations will be '
                        'applied. Use it at your own risk.\n'
                    )
                    break
        
        regex_pattern = regex.compile(user_regex)
        return(regex_pattern)

    elif legacy:
        # Keep the minimum value between `max_poly_a` provided and the remaining
        # length of the read after removing the barcode length.
        polya_run = min(max_poly_a, read2_length - longest_ab_tag - 1)

        # Read comment below for `s` meaning in the regex.
        pattern = r'(^(\L<options>)[TGC][A]{{{},}}){{{}<={}}}'.format(
            polya_run, error_type, maximum_distance
        )

    else:
        # `e` is part of the `regex` fuzzy logic: it means `error` in general,
        # whether it's a (s)ubstitution, (i)nsertion or (d)eletion. In this 
        # case, it means it allows `maximum_distance` errors to happen.
        pattern = r'(^\L<options>){{{}<={}}}'.format(error_type, maximum_distance)

    # Compiling the regex makes it run faster.
    regex_pattern = regex.compile(pattern, options=tag_keys)
    return(regex_pattern)

