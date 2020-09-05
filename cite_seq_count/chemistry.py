"""This module is holding code for all the remote fetching on the chemistries database."""
import requests
import sys
import os
import gzip
import io
import csv

from collections import namedtuple

from dataclasses import dataclass

from cite_seq_count import preprocessing

GLOBAL_LINK_RAW = "https://raw.githubusercontent.com/Hoohm/scg_lib_structs/10xv3_totalseq_b/chemistries/"
GLOBAL_LINK_GITHUB = "https://github.com/Hoohm/scg_lib_structs/raw/10xv3_totalseq_b/"
GLOBAL_LINK_GITHUB_IO = "https://teichlab.github.io/scg_lib_structs"
CHEMISTRY_DEFINITIONS = os.path.join(GLOBAL_LINK_RAW, "definitions.json")


@dataclass
class Chemistry:
    name: str
    cell_barcode_start: int
    cell_barcode_end: int
    umi_barcode_start: int
    umi_barcode_end: int
    R2_trim_start: int
    whitelist_path: str
    mapping_required: bool


def list_chemistries(url=CHEMISTRY_DEFINITIONS):
    """
    List all the available chemistries in the database
    Args:
        url (str): The url to the database file

    """
    print("Loading remote file from: {}".format(url))
    with requests.get(url) as r:
        r.raise_for_status()
    all_chemistry_defs = r.json()
    print(
        "Here are all the possible chemistries available at {}".format(
            GLOBAL_LINK_GITHUB_IO
        )
    )
    for chemistry in all_chemistry_defs:
        print(
            "\n-- {}\n   shortname: {}\n   Protocol link: {}\n\n".format(
                all_chemistry_defs[chemistry]["Description"],
                chemistry,
                os.path.join(
                    GLOBAL_LINK_GITHUB_IO, all_chemistry_defs[chemistry]["html"]
                ),
            )
        )


def get_chemistry_definition(chemistry_short_name, url=CHEMISTRY_DEFINITIONS):
    """
    Fetches chemistry definitions from a remote definitions.json and returns the json.
    """
    print("Loading remote file from: {}".format(url))
    with requests.get(url) as r:
        r.raise_for_status()
    chemistry_defs = r.json().get(chemistry_short_name, False)
    if not chemistry_defs:
        sys.exit(
            "Could not find the chemistry: {}. Please check that it does exist at: {}\nExiting".format(
                chemistry_short_name, url
            )
        )
    chemistry_def = Chemistry(
        name=chemistry_short_name,
        cell_barcode_start=chemistry_defs["barcode_structure_indexes"]["cell_barcode"][
            "R1"
        ]["start"],
        cell_barcode_end=chemistry_defs["barcode_structure_indexes"]["cell_barcode"][
            "R1"
        ]["stop"],
        umi_barcode_start=chemistry_defs["barcode_structure_indexes"]["umi_barcode"][
            "R1"
        ]["start"],
        umi_barcode_end=chemistry_defs["barcode_structure_indexes"]["umi_barcode"][
            "R1"
        ]["stop"],
        R2_trim_start=chemistry_defs["sequence_structure_indexes"]["R2"]["start"] - 1,
        whitelist_path=os.path.join(
            GLOBAL_LINK_GITHUB, "chemistries", chemistry_defs["whitelist"]["path"]
        ),
        mapping_required=chemistry_defs["whitelist"]["mapping"],
    )
    return chemistry_def


def get_csv_reader(file):
    if file.startswith("http://") or file.startswith("https://"):
        response = requests.get(file)
        response.raise_for_status()

        if file.endswith(".gz"):
            content = response.content
            text = gzip.decompress(content).decode("utf-8")
        else:
            text = response.text
        reader = csv.reader(io.StringIO(text))

    elif file.endswith(".gz"):
        f = gzip.open(file, mode="rt")
        reader = csv.reader(f)
    else:
        f = open(file, encoding="UTF-8")
        reader = csv.reader(f)

    return reader


def create_chemistry_definition(args):
    chemistry_def = Chemistry(
        name="custom",
        cell_barcode_start=args.cb_first,
        cell_barcode_end=args.cb_last,
        umi_barcode_start=args.umi_first,
        umi_barcode_end=args.umi_last,
        R2_trim_start=args.start_trim,
        whitelist_path=args.whitelist,
        mapping_required=args.translation,
    )
    return chemistry_def


def setup_chemistry(args):
    if args.chemistry:
        chemistry_def = get_chemistry_definition(args.chemistry)
        (whitelist, args.bc_threshold) = preprocessing.parse_whitelist_csv(
            csv_reader=get_csv_reader(chemistry_def.whitelist_path),
            barcode_length=chemistry_def.cell_barcode_end
            - chemistry_def.cell_barcode_start
            + 1,
            collapsing_threshold=args.bc_threshold,
        )
    else:
        chemistry_def = create_chemistry_definition(args)
        if args.whitelist:
            print("Loading whitelist")
            (whitelist, args.bc_threshold) = preprocessing.parse_whitelist_csv(
                csv_reader=get_csv_reader(args.whitelist),
                barcode_length=args.cb_last - args.cb_first + 1,
                collapsing_threshold=args.bc_threshold,
            )
        else:
            whitelist = False
    return (whitelist, chemistry_def)
