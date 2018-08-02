#!/usr/bin/env python3
"""
Authors: Christoph Hafemeister, Patrick Roelli
"""
import sys
import gzip
import csv
from collections import defaultdict
from collections import OrderedDict
from itertools import islice
from itertools import combinations
import pandas as pd
import time
import locale
import Levenshtein
import regex
import argparse
from argparse import RawTextHelpFormatter


def get_args():
    """
    Get args.
    """
    desc = "This script counts matching antobody tags from two fastq files."
    parser = argparse.ArgumentParser(prog='CITE Seq Count', description=desc,
                                     formatter_class=RawTextHelpFormatter)

    input_desc = "Required input files."
    inputs = parser.add_argument_group('Inputs', description=input_desc)

    inputs.add_argument('-R1', '--read1', dest='read1_path', required=True,
                        help="The path of read1 in gz format.")
    inputs.add_argument('-R2', '--read2', dest='read2_path', required=True,
                        help="The path of read2 in gz format.")
    inputs.add_argument('-t', '--tags', dest='tags', required=True,
                        help=("The path to the csv file containing the "
                              "antibody\nbarcodes as well as their respective "
                              "names.\n\nExample of an antibody barcode file "
                              "structure:\n\n"
                              "ATGCGA,First_tag_name\n"
                              "GTCATG,Second_tag_name"))

    bc_desc = ("Positions of the cellular barcodes and UMI.  If your cellular "
               "barcodes and UMI are positioned as follows:\n"
               "\tBarcodes from 1 to 16 and UMI from 17 to 26\n"
               "then this is the input you need:\n"
               "\t-cbf 1 -cbl 16 -umif 17 -umil 26")
    barcodes = parser.add_argument_group('Barcodes', description=bc_desc)

    barcodes.add_argument('-cbf', '--cell_barcode_first_base', dest='cb_first',
                          required=True, type=int,
                          help=("Postion of the first base of your cell "
                                "barcodes."))
    barcodes.add_argument('-cbl', '--cell_barcode_last_base', dest='cb_last',
                          required=True, type=int,
                          help=("Postion of the last base of your cell "
                                "barcodes."))
    barcodes.add_argument('-umif', '--umi_first_base', dest='umi_first',
                          required=True, type=int,
                          help="Postion of the first base of your UMI.")
    barcodes.add_argument('-umil', '--umi_last_base', dest='umi_last',
                          required=True, type=int,
                          help="Postion of the last base of your UMI.")

    barcodes_filtering = parser.add_mutually_exclusive_group(required=True)
    barcodes_filtering.add_argument('-cells', '--expected_cells', dest='cells',
                                    required=False, type=int,
                                    help=("Number of expected cells from your "
                                          "run."))
    whitelist_help = ("A csv file containning a whitelist of barcodes produced"
                      " by the mRNA data.\n\n\tExample:\n\tATGCTAGTGCTA\n"
                      "\tGCTAGTCAGGAT\n\tCGACTGCTAACG\n")
    barcodes_filtering.add_argument('-wl', '--whitelist', dest='whitelist',
                                    required=False, type=str,
                                    help=whitelist_help)

    filters_desc = ("Filtering for structure of antibody barcodes as well as "
                    "maximum hamming distance.")
    filters = parser.add_argument_group('filters', description=filters_desc)
    filters.add_argument('-tr', '--TAG_regex',
                        help="Only use if you know what you are doing."
                        "The regex that will be used to validate\n"
                        "an antibody barcode structure. Must be given in regex syntax."
                        "example:"
                        "\"^[ATGC]{6}[TGC][A]{6,}\"",
                        dest='tag_regex',
                        required=False,
                        type=str)
    filters.add_argument('-hd', '--hamming-distance', dest='hamming_thresh',
                         required=True, type=int,
                         help=("Maximum hamming distance allowed for antibody "
                               "barcode."))
    parser.add_argument('-n', '--first_n', required=False, type=int,
                        dest='first_n', default=None,
                        help="Select n reads to run on instead of all.")
    parser.add_argument('-o', '--output', required=True, type=str,
                        dest='outfile', help="Write result to file.")
    parser.add_argument('--debug', action='store_true',
                        help="Print extra information for debugging.")

    return parser


def parse_whitelist_csv(args):
    file = open(args.whitelist, mode='r')
    csvReader = csv.reader(file)
    length_barcodes = args.cb_last - args.cb_first + 1
    whitelist = [row[0].strip() for row in csvReader
                 if (len(row[0].strip()) == length_barcodes)]
    return set(whitelist)


def parse_tags_csv(filename):
    file = open(filename, mode='r')
    csvReader = csv.reader(file)
    odict = OrderedDict()
    for row in csvReader:
        odict[row[0]] = row[1]
    return odict

def check_tags(ab_map, maximum_dist):
    ab_barcodes = ab_map.keys()
    print(ab_barcodes)
    for a,b in combinations(ab_barcodes,2):
        if(Levenshtein.distance(a,b)<= maximum_dist):
            sys.exit('Minimum hamming distance of TAGS barcode is less than given threshold\nPlease use a smaller distance; exiting')


def generate_regex(ab_map, args, num_polyA):
    """Generate regex based ont he provided TAGS"""
    TAGS = ab_map.keys()
    lengths = OrderedDict()
    for TAG in TAGS:
        if (len(TAG) in lengths.keys()):
            lengths[len(TAG)]['mapping'][TAG]=ab_map[TAG]
        else:
            lengths[len(TAG)]=OrderedDict()
            lengths[len(TAG)]['mapping'] = OrderedDict()
            lengths[len(TAG)]['mapping'][TAG] = ab_map[TAG]
    #If there is only one length and the user provides a regex, us the users regex
    if ((len(lengths)==1) & (args.tag_regex is not None)):
        for length in lengths.keys():
            lengths[length]['regex'] = args.tag_regex
        return(lengths)
    if((len(lengths) != 1) & (args.tag_regex is not None)):
        exit('You cannot use your own regex with tag barcodes of different lengths')
    for length in lengths.keys():
        pattern = [''] * length
        for TAG in lengths[length]['mapping'].keys():
            for position in range(0,length):
                if (TAG[position] in pattern[position]):
                    continue
                else:
                    pattern[position] += TAG[position]
        lengths[length]['regex'] = '^([{}])[A]{{{},}}'.format(']['.join(pattern), num_polyA)
    return(lengths)



def main():
    parser = get_args()
    if not sys.argv[1:]:
        parser.print_help(file=sys.stderr)
        sys.exit(2)

    # Load args
    args = parser.parse_args()
    if args.whitelist:
        whitelist = parse_whitelist_csv(args)

    # Load TAGS barcodes
    ab_map = parse_tags_csv(args.tags)
    check_tags(ab_map, args.hamming_thresh)
    #Generate regex patterns auto
    regex_patterns = generate_regex(ab_map, args, num_polyA=6)

    
    # Create a set for UMI reduction. Fast way to check if it already exists
    UMI_reduce = set()
    # Create result table
    res_table = defaultdict(lambda: defaultdict(int))

    # Set counter
    n = 0
    top_n = None
    if args.first_n:
        top_n = args.first_n * 4

    # Define slices
    barcode_length = args.cb_last - args.cb_first + 1
    umi_length = args.umi_last - args.umi_first + 1
    barcode_umi_length = barcode_length + umi_length
    barcode_slice = slice(args.cb_first - 1, args.cb_last)
    umi_slice = slice(args.umi_first - 1, args.umi_last)

    unique_lines = set()
    with gzip.open(args.read1_path, 'rt') as textfile1, \
            gzip.open(args.read2_path, 'rt') as textfile2:
        # Read all 2nd lines from 4 line chunks
        secondlines = islice(zip(textfile1, textfile2), 1, top_n, 4)
        print('loading')

        t = time.time()
        for x, y in secondlines:
            x = x.strip()
            y = y.strip()
            line = x[barcode_slice] + x[umi_slice] + y
            unique_lines.add(line)

            n += 1
            if n % 1000000 == 0:
                print("Loaded last 1,000,000 lines in {:.3} seconds. Total "
                      "lines loaded {:,} ".format(time.time()-t, n))
                t = time.time()

        print('{} lines loaded'.format(n))
        print('{:,} uniques lines loaded'.format(len(unique_lines)))

        n = 0
        for line in unique_lines:
            cell_barcode = line[0:barcode_length]
            if args.whitelist:
                if cell_barcode not in whitelist:
                    continue

            UMI = line[barcode_length:barcode_umi_length]
            TAG_seq = line[barcode_umi_length:]
            BC_UMI_TAG = cell_barcode + UMI + TAG_seq
            if args.debug:
                print("{0}\t{1}\t{2}\t{3}".format(line, cell_barcode,
                                                  UMI, TAG_seq))

            # Check if UMI + TAG already in the set
            if BC_UMI_TAG not in UMI_reduce:
                # Check structure of the TAG
                no_structure_match=True
                for length in regex_patterns.keys():
                    match = regex.search(r'(?:({})){{i<={}}}'.format(regex_patterns[length]['regex'],args.hamming_thresh), TAG_seq)
                    if args.debug:
                        print("{0}\t{1}".format(regex_patterns[length]['regex'], TAG_seq))
                        print(match)

                    if match:
                        no_structure_match=False
                        TAG_seq = match.group(0)[0:length]

                        # Increment read count
                        res_table[cell_barcode]['total_reads'] += 1

                        # Get distance from all barcodes
                        temp_res = defaultdict()
                        for key, value in regex_patterns[length]['mapping'].items():
                            temp_res[value] = Levenshtein.hamming(TAG_seq, key)

                        # Get smallest value and get respective tag_name
                        min_value = min(temp_res.values())
                        min_index = list(temp_res.values()).index(min_value)
                        best = list(temp_res.keys())[min_index]

                        # ambiguous
                        if not isinstance(min_value, int):
                            if args.debug:
                                print("{0}\t{1}\t{2}\t{3}\t{4}".format(
                                      cell_barcode, UMI, x, y, TAG_seq))

                            res_table[cell_barcode]['ambiguous'] += 1
                            continue

                        # If over threshold
                        if min_value >= args.hamming_thresh:
                            res_table[cell_barcode]['no_match'] += 1
                            continue

                        res_table[cell_barcode][best] += 1

                        # Increment bad structure
                if(no_structure_match):
                    res_table[cell_barcode]['bad_struct'] += 1

                # Add BC_UMI_TAG to set
                UMI_reduce.add(BC_UMI_TAG)

            n += 1
            if n % 1000000 == 0:
                print("Processed 1,000,000 lines in {:.4} secondes. Total "
                      "lines processed: {:,}".format(time.time()-t, n))
                t = time.time()

    print("Done counting")
    res_matrix = pd.DataFrame(res_table)
    if ('total_reads' not in res_matrix.index):
        exit('No match found. Please check your regex or tags file')
    #Add potential missing cells if whitelist is used
    if(args.whitelist):
        res_matrix = res_matrix.reindex(whitelist, axis=1,fill_value=0)
    res_matrix.fillna(0, inplace=True)
    if args.cells:
        most_reads_ordered = res_matrix.sort_values(by='total_reads',
                                                    ascending=False,
                                                    axis=1).columns
        n_top_cells = int(args.cells + args.cells/100 * 30)
        top_Cells = most_reads_ordered[0:n_top_cells]
        res_matrix = res_matrix.loc[:, res_matrix.columns.isin(top_Cells)]

    res_matrix.to_csv(args.outfile, float_format='%.f')

if __name__ == '__main__':
    main()
