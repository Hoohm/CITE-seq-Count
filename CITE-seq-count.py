
#! /usr/bin/python3
#Authors: Christoph Hafemeister, Patrick Roelli
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
import distance
import re
import argparse
from argparse import RawTextHelpFormatter

def get_args():
    """Get args."""
    parser = argparse.ArgumentParser(prog='CITE Seq Count',
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
                        help="""The path to the csv file containing the TAGS\nbarcodes as well as their respective names.

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
    barcodes.add_argument('-cells', '--expected_cells',
                        help='Number of expected cells from your run',
                        dest='cells',
                        required=False, 
                        type=int)
    
    filters = parser.add_argument_group('filters', description="""Filtering for structure of TAGS as well as maximum hamming distance.""")
    filters.add_argument('-tr', '--TAG_regex',
                        help="""The regex that will be used to validate a tag structure. Must be given in regex syntax.
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
                        required=True,
                        help='write to file instead of stdout',
                        type=str,
                        dest='outfile')
    return parser


def parse_csv(filename):
    file = open(filename, mode='r')
    csvReader = csv.reader(file)
    odict = OrderedDict()
    for row in csvReader:
        #odict[row[0].encode('utf-8')] = row[1]
        odict[row[0]] = row[1]
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
    #TAG_structure = re.compile(args.tag_regex.encode('utf-8'))
    TAG_structure = re.compile(args.tag_regex)
    #Create a set for UMI reduction. Fast way to check if it already exists
    UMI_reduce = set()
    #Create result table
    res_table = defaultdict(lambda : defaultdict(int))
    # set counter
    n = 0
    unique_lines = set()
    with gzip.open(args.read1_path, 'rt') as textfile1, gzip.open(args.read2_path, 'rt') as textfile2: 
        #Read all 2nd lines from 4 line chunks
        secondlines = islice(zip(textfile1, textfile2), 1, args.first_n, 4)
        print('loading')
        t = time.time()
        for x, y in secondlines:
            #print(x,y)
            x = x.strip()
            y = y.strip()
            #print(x,y)
            line = x[(args.cb_first-1):args.cb_last] +x[args.umi_first-1: args.umi_last] + y
            unique_lines.add(line)
            n+=1
            if(n%1000000 == 0):
                print("Loaded last 1,000,000 lines in {:.3} seconds. Total lines loaded {:,} ".format(time.time()-t, n))
                t = time.time()
        print('{} lines loaded'.format(n))
        print('{:,} uniques lines loaded'.format(len(unique_lines)))
        n=0
        for line in unique_lines:
            cell_barcode = line[0:args.cb_last-args.cb_first+1]
            UMI = line[len(cell_barcode):len(cell_barcode)+args.umi_last-args.umi_first+1]
            BC_UMI = cell_barcode + UMI
            TAG_seq = line[len(BC_UMI):]
            BC_UMI_TAG = cell_barcode + UMI + TAG_seq
            #print("{0}\t{1}\t{2}\t{3}".format(line, cell_barcode, UMI,TAG_seq))
            if BC_UMI_TAG not in UMI_reduce:#check if UMI + TAG already in the set
                if re.match(TAG_structure, TAG_seq):#check structure of the TAG
                    TAG_seq = TAG_seq[0:tag_length]
                    res_table[cell_barcode]['total_reads'] += 1 #increment read count
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
                    res_table[cell_barcode][best]+=1
                else:
                    res_table[cell_barcode]['bad_struct'] += 1 #Increment bad structure
                UMI_reduce.add(BC_UMI_TAG) # Add BC_UMI_TAG to set    
            n+=1
            if(n%1000000 == 0):
                print("Processed 1,000,000 lines in {:.4} secondes. Total lines processed: {:,}".format(time.time()-t, n))
                t = time.time()
    print('Done counting')
    
    res_matrix = pd.DataFrame(res_table)
    res_matrix.fillna(0, inplace=True)
    res_matrix.rename(columns={'': 'Tags'}, inplace=True)
    print(res_matrix.columns)

    most_reads_ordered = res_matrix.sort_values(by='total_reads', ascending=False, axis=1).axes[1]
    n_top_cells = int(args.cells + args.cells/100 * 30)
    top_Cells = most_reads_ordered[0:(n_top_cells)]
    
    print(top_Cells)
    res_matrix.loc[:,(res_matrix.columns.isin(top_Cells))].to_csv(args.outfile, float_format='%.f')
if __name__ == '__main__':
    main()