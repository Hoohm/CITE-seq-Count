
#! /usr/binp/ython3
import sys
import gzip
import csv
from collections import defaultdict
from collections import OrderedDict
from itertools import islice
from itertools import combinations

import distance
import re
import argparse
from argparse import RawTextHelpFormatter

def get_args():
    """Get args."""
    parser = argparse.ArgumentParser(prog='count_tags',
        description="""This script counts matching tags from two fastq files.""", formatter_class=RawTextHelpFormatter)
    inputs = parser.add_argument_group('Inputs')
    inputs.add_argument('-R1', '--read1',
                        action='store',
                        help='The path of read1 in gz format',
                        dest='read1_path',
                        required=True)
    inputs.add_argument('-R2', '--read2',
                        help='The path of read2 in gz format',
                        dest='read2_path',
                        required=True)
    inputs.add_argument('-t', '--tags',
                        help="""The path to the csv file containing the TAGS\narcodes as well as their respective names.

Example of TAG file structure:

ATGCGA,First_tag_name
GTCATG,Second_tag_name""",
                        dest='tags',
                        required=True)
    barcodes = parser.add_argument_group('Barcodes', description="""Positions of the cellular barcodes and UMI.
If your cellular barcodes and UMI are positioned as follows:
Barcodes from 1 to 16 and UMI from 17 to 26, then this is the input you need:
-cbf 1 -cbl 16 -umif 17 -umil 26""")
    barcodes.add_argument('-cbf', '--cell_barcode_first_base',
                        help='Postion of the first base of your cell barcodes',
                        dest='cb_first',
                        required=True, 
                        type=int)
    barcodes.add_argument('-cbl', '--cell_barcode_last_base',
                        help='Postion of the last base of your cell barcodes',
                        dest='cb_last',
                        required=True, 
                        type=int)
    barcodes.add_argument('-umif', '--umi_first_base',
                        help='Postion of the first base of your UMI',
                        dest='umi_first',
                        required=True, 
                        type=int)
    barcodes.add_argument('-umil', '--umi_last_base',
                        help='Postion of the last base of your UMI',
                        dest='umi_last',
                        required=True, 
                        type=int)
    filters = parser.add_argument_group('filters', description="""Filtering for structure of TAGS as well as maximum hamming distance.""")
    filters.add_argument('-tr', '--TAG_regex',
                        help="""The regex that will be sued to validate a tag structure. Must be given in regex syntax.
example:
\"^[ATGC]{6}[TGC][A]{6,}\"""",
                        dest='tag_regex',
                        required=True, 
                        type=str)
    filters.add_argument('-hd', '--hamming-distance',
                        help='Maximum hamming distance allowed',
                        dest='hamming_thresh',
                        required=True, 
                        type=int)
    parser.add_argument('-n','--first_n',
                        required=False,
                        help='Select n lines to run on instead of all.\nIf you want the first 10 reads, use -n 40',
                        type=int,
                        dest='first_n',
                        default=None)
    parser.add_argument('-o','--output',
                        required=False,
                        help='write to file instead of stdout',
                        type=str,
                        dest='outfile')
    return parser


def parse_csv(filename):
    file = open(filename, mode='r')
    csvReader = csv.reader(file)
    odict = OrderedDict()
    for row in csvReader:
        odict[row[0].encode('utf-8')] = row[1]
    return odict


def check_tags(ab_map, maximum_dist):
    ab_barcodes = ab_map.keys()
    for a,b in combinations(ab_barcodes,2):
        if(len(a)!=len(b)):
            sys.exit('Lenght of {} is different than length of {}. Can only run with all TAGS having the same length; exiting'.format(ab_map[a], ab_map[b]))
        if(distance.hamming(a,b)<= maximum_dist):
            sys.exit('Minimum hamming distance of TAGS barcode is less than given threshold\nPlease use a smaller distance; exiting')
    #Return length of TAGS. Since all are the same lenght, return the length of last a
    return(len(a))

def main():
    parser = get_args()
    if not sys.argv[1:]:
        sys.exit(parser.print_help())
    #Load args
    args = parser.parse_args()
    #Load TAGS barcodes
    ab_map = parse_csv(args.tags)
    #Chekck hamming threshold
    tag_length = check_tags(ab_map, args.hamming_thresh)
    #Create TAG structure filter
    TAG_structure = re.compile(args.tag_regex.encode('utf-8'))
    #Create a set for UMI reduction. Fast way to check if it already exists
    UMI_reduce = set()
    #Create result table
    res_table = defaultdict(lambda : defaultdict(int))
    with gzip.open(args.read1_path, 'r') as textfile1, gzip.open(args.read2_path, 'r') as textfile2: 
        #Read all 2nd lines from 4 line chunks
        secondlines = islice(zip(textfile1, textfile2), 1, args.first_n, 4)
        for (x, y) in secondlines:
            x = x.strip()
            y = y.strip()
            cell_barcode = x[(args.cb_first-1):args.cb_last]
            UMI = x[(args.umi_first-1):args.umi_last]
            BC_UMI = cell_barcode + UMI
            BC_UMI_TAG = BC_UMI + y[0:tag_length]
            TAG_seq = y[0:tag_length]
            #print("{0}\t{1}\t{2}\t{3}\t{4}".format(cell_barcode, UMI, x,y,TAG_seq))
            #Check if cell barcode has been found before. If not, create result structure
            if cell_barcode not in res_table:
                for TAG_name in ab_map.values():
                    res_table[cell_barcode][TAG_name]
                res_table[cell_barcode]['no_match']
                res_table[cell_barcode]['total_reads']
                res_table[cell_barcode]['bad_struct']
                res_table[cell_barcode]['ambiguous']
            if re.match(TAG_structure, y):#check structure of the TAG
                if BC_UMI_TAG in UMI_reduce:#check if UMI + TAG already in the set
                    continue#go to next y read
                else:
                    res_table[cell_barcode]['total_reads'] += 1 #increment read count
                    UMI_reduce.add(BC_UMI_TAG) # Add to set
                temp_res = defaultdict()
                for key, value in ab_map.items():
                    temp_res[value] = distance.hamming(TAG_seq, key) #Get distance from all barcodes
                best = list(temp_res.keys())[list(temp_res.values()).index(min(temp_res.values()))]#Get smallest value and get respective tag_name
                if(not isinstance(min(temp_res.values()),int)):# ambiguous
                    #print("{0}\t{1}\t{2}\t{3}\t{4}".format(cell_barcode, UMI, x,y,TAG_seq))
                    res_table[cell_barcode]['ambiguous'] += 1
                    continue #next entry
                if(min(temp_res.values()) >= args.hamming_thresh):#If over threshold
                    res_table[cell_barcode]['no_match'] += 1
                    continue #next entry
                res_table[cell_barcode][best] += 1
            else:
                res_table[cell_barcode]['bad_struct'] += 1 #Increment bad structure
    # Create header
    out_str = "{}\t{}\t{}\t{}\t{}\t{}\n".format('cell', "\t".join([x for x in ab_map.values()]), 'no_match', 'ambiguous','total_reads','bad_struct')
    # fill up result string
    for cell_barcode in res_table:
      out_str+=("{}\t{}\t{}\t{}\t{}\t{}\n".format(
        cell_barcode.decode("utf-8"),
        "\t".join([str(res_table[cell_barcode][x]) for x in ab_map.values()]),
        res_table[cell_barcode]['no_match'],
        res_table[cell_barcode]['ambiguous'],
        res_table[cell_barcode]['total_reads'],
        res_table[cell_barcode]['bad_struct']))
    if(args.outfile):
        with open(args.outfile, 'w') as h:
            h.write(out_str)
    else:
        sys.stdout.write(out_str)

if __name__ == '__main__':
    main()