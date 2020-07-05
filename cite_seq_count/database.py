"""This module is holding code for all the remote fetching on the chemistries database."""
import requests
import sys

from collections import namedtuple


CHEMISTRY_DEFINITIONS = "https://raw.githubusercontent.com/Hoohm/scg_lib_structs/10xv3_totalseq_b/chemistries/definitions.json"

CHEMISTRY_DEFINITION = namedtuple('chemistry_def', ['barcode_start','barcode_end','umi_start','umi_end', 'R2_trim_start'])

def get_chemistry_definition(chemistry_short_name, url=CHEMISTRY_DEFINITIONS):
    """
    Fetches chemistry definitions from a remote definitions.json and returns the json.
    """
    print('Loading remote file from: {}'.format(url))
    with requests.get(url) as r:
        r.raise_for_status()
    chemistry_defs = r.json().get(chemistry_short_name, False)
    if not chemistry_defs:
        sys.exit('Could not find the chemistry: {}. Please check that it does exist at: {}\nExiting'.format(chemistry_short_name, url))
    
    chemistry_def = CHEMISTRY_DEFINITION(
        barcode_start=chemistry_defs['barcode_structure_indexes']['cell_barcode']['R1']['start'],
        barcode_end=chemistry_defs['barcode_structure_indexes']['cell_barcode']['R1']['stop'],
        umi_start=chemistry_defs['barcode_structure_indexes']['umi_barcode']['R1']['start'],
        umi_end=chemistry_defs['barcode_structure_indexes']['umi_barcode']['R1']['stop'],
        R2_trim_start=chemistry_defs['sequence_structure_indexes']['R2']['start'] - 1
    )
    return(chemistry_def)

def create_chemistry_definition(args):
    """
    """
    chemistry_def = CHEMISTRY_DEFINITION(
        barcode_start=args.cbf,
        barcode_end=args.cbl,
        umi_start=args.umif,
        umi_end=args.umil,
        R2_trim_start=args.start_trim
    )
    return(chemistry_def)