import csv
import gzip
import sys
import regex
import Levenshtein

from math import floor
from collections import OrderedDict
from itertools import combinations
from itertools import islice 

def get_indexes(start_index, chunk_size, nth):
    """
    Creates indexes from a reference index, a chunk size an nth number

    Args:
        start_index (int): first position
        chunk_size (int): Chunk size
        nth (int): The nth number
    
    Returns:
        list: First and last position of indexes
    """
    start_index = nth * chunk_size
    stop_index = chunk_size + nth * chunk_size
    return([start_index,stop_index])


def chunk_reads(n_reads, n):
    """
    Creates a list of indexes for the islice iterator from the map_reads function.

    Args:
        n_reads (int): Number of reads to split
        n (int): How many buckets for the split.
    Returns:
        indexes (list(list)): Each entry contains the first and the last index for a read.
    """
    indexes=list()
    if n_reads % n == 0:
        chunk_size = int(n_reads/n)
        rest = 0
    else:
        chunk_size = floor(n_reads/n)
        rest = n_reads - (n*chunk_size)
    for i in range(0,n):
        indexes.append(get_indexes(i, chunk_size, i))
    indexes[-1][1] += rest
    return(indexes)
    

def parse_whitelist_csv(filename, barcode_length, collapsing_threshold):
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
        cell_pattern = regex.compile(r'[ATGC]{{{}}}'.format(barcode_length))
        whitelist = [row[0].strip(STRIP_CHARS) for row in csv_reader
                     if (len(row[0].strip(STRIP_CHARS)) == barcode_length)]
    for cell_barcode in whitelist:
        if not cell_pattern.match(cell_barcode):
            sys.exit('This barcode {} is not only composed of ATGC bases.'.format(cell_barcode))
    collapsing_threshold=test_cell_distances(whitelist, collapsing_threshold)
    return(set(whitelist), collapsing_threshold)


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
        print('Testing cell barcode collapsing threshold of {}'.format(collapsing_threshold))
        for comb in combinations(whitelist, 2):
            if Levenshtein.hamming(comb[0], comb[1]) <= collapsing_threshold:
                collapsing_threshold -= 1
                print('Value is too high, reducing it by 1')
                break
        else:
            ok = True
    print('Using {} for cell barcode collapsing threshold'.format(collapsing_threshold))
    return(collapsing_threshold)


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
        if (distance <= (maximum_distance - 1)):
            offending_pairs.append([a, b, distance])
    DNA_pattern = regex.compile('^[ATGC]*$')
    for tag in tags:
        if not DNA_pattern.match(tag):
            print('This tag {} is not only composed of ATGC bases.\nPlease check your tags file'.format(tag))
            sys.exit('Exiting the application.\n')
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
    print('Counting number of reads')
    with gzip.open(file_path, "rt",encoding="utf-8",errors='ignore') as f:
        n_lines = sum(bl.count("\n") for bl in blocks(f))
    if n_lines %4 !=0:
        sys.exit('{}\'s number of lines is not a multiple of 4. The file might be corrupted.\n Exiting')
    return(n_lines)

