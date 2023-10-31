"""This module is holding code for all the remote fetching on the chemistries database."""
import os
import pooch
import json


from dataclasses import dataclass
from argparse import ArgumentParser
from cite_seq_count import preprocessing
import polars as pl

GLOBAL_LINK_RAW = "https://raw.githubusercontent.com/Hoohm/scg_lib_structs/10xv3_totalseq_b/chemistries/"
GLOBAL_LINK_GITHUB = "https://github.com/Hoohm/scg_lib_structs/raw/10xv3_totalseq_b/"
GLOBAL_LINK_GITHUB_IO = "https://teichlab.github.io/scg_lib_structs"
# "https://github.com/Hoohm/scg_lib_structs/raw/10xv3_totalseq_b/chemistries/whitelists/3M-february-2018.csv.gz"
# CHEMISTRY_DEFINITIONS = os.path.join(GLOBAL_LINK_RAW, "definitions.json")


@dataclass
class Chemistry:
    name: str
    cell_barcode_start: int
    cell_barcode_end: int
    umi_barcode_start: int
    umi_barcode_end: int
    r2_trim_start: int
    barcode_reference_path: str


DEFINITIONS_DB = pooch.create(
    path=pooch.os_cache("cite_seq_count"),
    base_url=GLOBAL_LINK_RAW,
    version="0.1.0",
    env="MYPACKAGE_DATA_DIR",
    # The cache file registry. A dictionary with all files managed by this
    # pooch. Keys are the file names (relative to *base_url*) and values
    # are their respective SHA256 hashes. Files will be downloaded
    # automatically when needed (see fetch_gravity_data).
    registry={
        "definitions.json": "4f2e1cee60446062f0c805b0b51ce12318cac04c64f1c7825037e2c73ff9955e"
    },
)


def fetch_definitions() -> dict:
    """
    Load some sample gravity data to use in your docs.
    """
    fname = DEFINITIONS_DB.fetch("definitions.json")
    with open(fname, encoding="utf-8") as json_file:
        data = json_file.read()
    json_data = json.loads(data)
    return json_data


def list_chemistries(all_chemistry_defs: str) -> None:
    """
    List all the available chemistries in the database
    Args:
        url (str): The url to the database file

    """
    all_chemistry_defs = fetch_definitions()
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


def get_chemistry_definition(chemistry_short_name: str) -> Chemistry:
    """
    Fetches chemistry definitions from a remote definitions.json and returns the json.
    """
    chemistry_defs = fetch_definitions()[chemistry_short_name]
    print(chemistry_defs)
    if chemistry_defs["barcode_reference"]["path"] not in DEFINITIONS_DB.registry:
        path = pooch.retrieve(
            url=os.path.join(
                GLOBAL_LINK_GITHUB,
                "chemistries",
                chemistry_defs["barcode_reference"]["path"],
            ),
            known_hash=None,
            fname=chemistry_defs["barcode_reference"]["path"],
            path=DEFINITIONS_DB.abspath,
        )
    else:
        path = DEFINITIONS_DB.registry[chemistry_defs["barcode_reference"]["path"]]
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
        r2_trim_start=chemistry_defs["sequence_structure_indexes"]["R2"]["start"] - 1,
        barcode_reference_path=path,
    )
    return chemistry_def


def create_chemistry_definition(args: ArgumentParser) -> Chemistry:
    chemistry_def = Chemistry(
        name="custom",
        cell_barcode_start=args.cb_first,
        cell_barcode_end=args.cb_last,
        umi_barcode_start=args.umi_first,
        umi_barcode_end=args.umi_last,
        r2_trim_start=args.start_trim,
        barcode_reference_path=args.barcode_reference,
    )
    return chemistry_def


def setup_chemistry(args: ArgumentParser) -> tuple[pl.DataFrame | None, Chemistry]:
    if args.chemistry_id:
        chemistry_def = get_chemistry_definition(args.chemistry_id)
        barcode_reference = preprocessing.parse_barcode_reference(
            filename=chemistry_def.barcode_reference_path,
            barcode_length=chemistry_def.cell_barcode_end
            - chemistry_def.cell_barcode_start
            + 1,
            required_header=["reference"],
        )
    else:
        chemistry_def = create_chemistry_definition(args)
        if args.barcode_reference:
            print("Loading barcode reference")
            barcode_reference = preprocessing.parse_barcode_reference(
                filename=args.barcode_reference,
                barcode_length=args.cb_last - args.cb_first + 1,
                required_header=["reference"],
            )
        else:
            barcode_reference = None
    return (barcode_reference, chemistry_def)
