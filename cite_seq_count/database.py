"""This module is holding code for all the remote fetching on the chemistries database."""
import requests
import sys
import os

from collections import namedtuple
from tempfile import TemporaryFile
from dataclasses import dataclass

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


def create_chemistry_definition(args):
    return chemistry_def
