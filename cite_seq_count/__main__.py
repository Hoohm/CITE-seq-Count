#!/usr/bin/env python3
"""
Authors: Christoph Hafemeister, Patrick Roelli
"""
import sys
import gzip
import csv
import warnings
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
import pkg_resources
version = pkg_resources.require("cite_seq_count")[0].version


def get_args():
    """
    Get args.
    """
    desc = "This script counts matching antobody tags from two fastq files. Version {}".format(version)
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
    filters.add_argument('-hd', '--hamming-distance', dest='hamming_thresh',
                         required=False, type=int, default=2,
                         help=("Maximum hamming distance allowed for antibody "
                               "barcode."))
    parser.add_argument('-n', '--first_n', required=False, type=int,
                        dest='first_n', default=None,
                        help="Select n reads to run on instead of all.")
    parser.add_argument('-o', '--output', required=True, type=str,
                        dest='outfile', help="Write result to file.")
    parser.add_argument('-u', '--unknown-tags', required=False, type=str,
                        dest='unknowns_file', help="Write table of unknown tags to file.")
    parser.add_argument('--debug', action='store_true',
                        help="Print extra information for debugging.")
    regex_pattern = parser.add_mutually_exclusive_group(required=False)
    regex_pattern.add_argument('-tr', '--TAG_regex',
                        help="Only use if you know what you are doing."
                        "The regex that will be used to validate\n"
                        "an antibody barcode structure. Must be given in regex syntax."
                        "example:"
                        "\"^[ATGC]{6}[TGC][A]{6,}\"",
                        dest='tag_regex',
                        required=False,
                        type=str)
    regex_pattern.add_argument('-l', '--legacy', required=False,
                        dest='legacy', default=False, action='store_true',
                        help="Use this option if you used an earlier versions"
                        " of the kit that adds a T C or G at the end and you"
                        " expect polyA tails in the data.")
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
        odict[row[0].strip()] = row[1].strip()
    return odict

def check_tags(ab_map, maximum_dist):
    # Adding the barcode to the name of the TAG
    # This means we don't need to share the mapping of the antibody and the barcode.
    new_ab_map = {}
    for TAG in ab_map:
        new_ab_map[TAG] = ab_map[TAG]  + '-' + TAG
    if(len(ab_map) == 1):
        return(new_ab_map)
    for a,b in combinations(new_ab_map.keys(),2):
        if(Levenshtein.distance(a,b)<= maximum_dist):
            sys.exit('Minimum hamming distance of TAGS barcode is less than given threshold\nPlease use a smaller distance; exiting')
    return(new_ab_map)

def generate_regex(ab_map, args, R2_length, max_polyA):
    """Generate regex based ont he provided TAGS"""
    lengths = OrderedDict()
    for TAG in ab_map:
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
        if(args.legacy):
            lengths[length]['regex'] = '^([{}])[TGC][A]{{{},}}'.format(']['.join(pattern), min(max_polyA,(R2_length-length-1)))
        else:
            lengths[length]['regex'] = '^([{}])'.format(']['.join(pattern))
    return(lengths)

def get_read_length(file_path):
    with gzip.open(file_path, 'r') as fastq_file:
        secondlines = islice(fastq_file, 1, 1000, 4)
        temp_length = len(next(secondlines).rstrip())
        for sequence in secondlines:
            read_length = len(sequence.rstrip())
            if(temp_length != read_length):
                sys.exit('Read2 length is not consistent, please trim all Read2 reads at the same length')
            temp_length = read_length
    return(read_length)


def check_read_lengths(R1_length, R2_length, args):
    barcode_length = args.cb_last - args.cb_first + 1
    umi_length = args.umi_last - args.umi_first + 1
    barcode_umi_length = barcode_length + umi_length
    barcode_slice = slice(args.cb_first - 1, args.cb_last)
    umi_slice = slice(args.umi_first - 1, args.umi_last)

    if(barcode_umi_length) > R1_length:
        sys.exit('Read 1 length is shorter than the option you are using for cell and UMI barcodes length. Please check your options and rerun.')
    elif(barcode_umi_length) < R1_length:
        print('**WARNING**\nRead 1 length is {}bp but you are using {}bp for cell and UMI barcodes combined.\nThis might lead to wrong cell attribution and skewed umi counts.\n'.format(R1_length, barcode_umi_length))
    return(barcode_slice, umi_slice, barcode_umi_length)

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
    ab_map = check_tags(ab_map, args.hamming_thresh)
    
    #Get read lengths
    R1_length = get_read_length(args.read1_path)
    R2_length = get_read_length(args.read2_path)

    #Generate regex patterns automatically
    regex_patterns = generate_regex(ab_map=ab_map, args=args, R2_length=R2_length, max_polyA=6)
    if(args.debug):
        print(regex_patterns)
    
    # Create a set for UMI reduction. Fast way to check if it already exists
    UMI_reduce = set()
    # Create result table
    res_table = defaultdict(lambda: defaultdict(int))
    no_match_table = defaultdict(int)

    # Set counter
    n = 0
        
    # Check that read 1 and options match and define slices
    (barcode_slice, umi_slice, barcode_umi_length) = check_read_lengths(R1_length, R2_length, args)
    
    unique_lines = set()
    with gzip.open(args.read1_path, 'rt') as textfile1, \
            gzip.open(args.read2_path, 'rt') as textfile2:
        # Read all 2nd lines from 4 line chunks. If first_n not None read only 4 times the given amount.
        secondlines = islice(zip(textfile1, textfile2), 1, (args.first_n * 4 if args.first_n is not None else args.first_n), 4)
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

        print('{:,} reads loaded'.format(n))
        print('{:,} uniques reads loaded'.format(len(unique_lines)))

        n = 1
        for line in unique_lines:
            if n % 1000000 == 0:
                print("Processed 1,000,000 lines in {:.4} secondes. Total "
                      "lines processed: {:,}".format(time.time()-t, n))
                t = time.time()

            cell_barcode = line[barcode_slice]
            if args.whitelist:
                if cell_barcode not in whitelist:
                    n += 1
                    continue


            UMI = line[umi_slice]
            TAG_seq = line[barcode_umi_length:]
            BC_UMI_TAG = cell_barcode + UMI + TAG_seq
            if args.debug:
                print("\nline:{0}\ncell_barcode:{1}\tUMI:{2}\tTAG_seq:{3}\nline length:{4}\tcell barcode length:{5}\tUMI length:{6}\tTAG sequence length:{7}".format(line, cell_barcode,
                                                  UMI, TAG_seq, len(line), len(cell_barcode), len(UMI), len(TAG_seq)))

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
                            no_match_table[TAG_seq] += 1
                            continue

                        res_table[cell_barcode][best] += 1

                        # Increment bad structure
                if(no_structure_match):
                    res_table[cell_barcode]['bad_struct'] += 1

                # Add BC_UMI_TAG to set
                UMI_reduce.add(BC_UMI_TAG)

            n += 1
            
    print("Done counting")
    res_matrix = pd.DataFrame(res_table)
    if ('total_reads' not in res_matrix.index):
        exit('No match found. Please check your regex or tags file')
    #Add potential missing cells if whitelist is used
    if args.whitelist:
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
    
    if args.unknowns_file:
        keys = list(no_match_table.keys())
        vals = list(no_match_table.values())
        no_match_matrix = pd.DataFrame({"tag": keys, "total": vals})
        no_match_matrix = no_match_matrix.sort_values(by='total', ascending=False)            
        no_match_matrix.to_csv(args.unknowns_file, float_format='%.f', index=False)
                  

if __name__ == '__main__':
    main()
