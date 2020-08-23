import pkg_resources
import sys


from argparse import ArgumentParser, ArgumentTypeError, RawTextHelpFormatter

# pylint: disable=no-name-in-module
from multiprocess import cpu_count


def get_package_version():
    version = pkg_resources.require("cite_seq_count")[0].version
    return version


def chunk_size_limit(arg):
    """Validates chunk_size limits"""
    max_size = 2147483647
    try:
        f = int(arg)
    except ValueError:
        raise ArgumentTypeError("Must be an int")
    if f < 1 or f > max_size:
        raise ArgumentTypeError(
            "Argument must be < " + str(max_size) + "and > " + str(1)
        )
    else:
        return False
    return f


def get_args():
    """
    Get args.
    """

    parser = ArgumentParser(
        prog="CITE-seq-Count",
        formatter_class=RawTextHelpFormatter,
        description=(
            "This script counts matching antibody tags from paired fastq "
            "files. Version {}".format(get_package_version())
        ),
    )

    # REQUIRED INPUTS group.
    inputs = parser.add_argument_group("Inputs", description="Required input files.")
    inputs.add_argument(
        "-R1",
        "--read1",
        dest="read1_path",
        required=True,
        help=(
            "The path of Read1 in gz format, or a comma-separated list of paths to all Read1 files in"
            " gz format (E.g. A1.fq.gz,B1.fq,gz,..."
        ),
    )
    inputs.add_argument(
        "-R2",
        "--read2",
        dest="read2_path",
        required=True,
        help=(
            "The path of Read2 in gz format, or a comma-separated list of paths to all Read2 files in"
            " gz format (E.g. A2.fq.gz,B2.fq,gz,..."
        ),
    )
    inputs.add_argument(
        "-t",
        "--tags",
        dest="tags",
        required=True,
        help=(
            "The path to the csv file containing the antibody\n"
            "barcodes as well as their respective names.\n\n"
            "Example of an antibody barcode file structure:\n\n"
            "\tATGCGA,First_tag_name\n"
            "\tGTCATG,Second_tag_name"
        ),
    )
    # BARCODES group.
    barcodes = parser.add_argument_group(
        "Barcodes",
        description=(
            "Positions of the cellular barcodes and UMI. If your "
            "cellular barcodes and UMI\n are positioned as follows:\n"
            "\tBarcodes from 1 to 16 and UMI from 17 to 26\n"
            "then this is the input you need:\n"
            "\t-cbf 1 -cbl 16 -umif 17 -umil 26"
        ),
    )
    barcodes.add_argument("--chemistry", type=str, required=False, default=False)
    if "--chemistry" not in sys.argv:
        barcodes.add_argument(
            "-cbf",
            "--cell_barcode_first_base",
            dest="cb_first",
            required=True,
            type=int,
            help=("Postion of the first base of your cell " "barcodes."),
        )
        barcodes.add_argument(
            "-cbl",
            "--cell_barcode_last_base",
            dest="cb_last",
            required=True,
            type=int,
            help=("Postion of the last base of your cell " "barcodes."),
        )
        barcodes.add_argument(
            "-umif",
            "--umi_first_base",
            dest="umi_first",
            required=True,
            type=int,
            help="Postion of the first base of your UMI.",
        )
        barcodes.add_argument(
            "-umil",
            "--umi_last_base",
            dest="umi_last",
            required=True,
            type=int,
            help="Postion of the last base of your UMI.",
        )
    barcodes.add_argument(
        "--umi_collapsing_dist",
        dest="umi_threshold",
        required=False,
        type=int,
        default=2,
        help="threshold for umi collapsing.",
    )
    barcodes.add_argument(
        "--no_umi_correction",
        required=False,
        action="store_true",
        default=False,
        dest="no_umi_correction",
        help="Deactivate UMI collapsing",
    )
    barcodes.add_argument(
        "--bc_collapsing_dist",
        dest="bc_threshold",
        required=False,
        type=int,
        default=1,
        help="threshold for cellular barcode collapsing.",
    )
    # Cells group
    cells = parser.add_argument_group(
        "Cells", description=("Expected number of cells and potential whitelist")
    )

    cells.add_argument(
        "-cells",
        "--expected_cells",
        dest="expected_cells",
        required=True,
        type=int,
        help=("Number of expected cells from your run."),
        default=0,
    )
    if "--chemistry" not in sys.argv:
        cells.add_argument(
            "-wl",
            "--whitelist",
            dest="whitelist",
            required=False,
            type=str,
            help=(
                "A csv file containning a whitelist of barcodes produced"
                " by the mRNA data.\n\n"
                "\tExample:\n"
                "\tATGCTAGTGCTA\n\tGCTAGTCAGGAT\n\tCGACTGCTAACG\n\n"
                "Or 10X-style:\n"
                "\tATGCTAGTGCTA-1\n\tGCTAGTCAGGAT-1\n\tCGACTGCTAACG-1\n"
            ),
        )

        cells.add_argument(
            "--translation",
            required=False,
            type=str,
            help="A csv file containing the mapping between two sets of cell barcode list.\n"
            "A required header such as the reference is named whitelist. Example:\n\n"
            "\twhitelist,trasnlation\n"
            "\tAAACCCAAGAAACACT,AAACCCATCAAACACT\n"
            "\tAAACCCAAGAAACCAT,AAACCCATCAAACCAT\n"
            "\nThe output matrix will possess both cell barcode IDs",
        )

    # FILTERS group.
    filters = parser.add_argument_group(
        "TAG filters", description=("Filtering and trimming for read2.")
    )
    filters.add_argument(
        "--max-errors",
        dest="max_error",
        required=False,
        type=int,
        default=2,
        help=("Maximum Levenshtein distance allowed for antibody barcodes."),
    )

    filters.add_argument(
        "-trim",
        "--start-trim",
        dest="start_trim",
        required=False,
        type=int,
        default=0,
        help=("Number of bases to discard from read2."),
    )

    filters.add_argument(
        "--sliding-window",
        dest="sliding_window",
        required=False,
        default=False,
        action="store_true",
        help=("Allow for a sliding window when aligning."),
    )

    # Parallel group.
    parallel = parser.add_argument_group(
        "Parallelization options",
        description=("Options for performance on parallelization"),
    )
    # Remaining arguments.
    parallel.add_argument(
        "-T",
        "--threads",
        required=False,
        type=int,
        dest="n_threads",
        default=cpu_count(),
        help="How many threads are to be used for running the program",
    )
    parallel.add_argument(
        "-C",
        "--chunk_size",
        required=False,
        type=chunk_size_limit,
        dest="chunk_size",
        help="How many reads should be sent to a child process at a time",
    )
    parallel.add_argument(
        "--temp_path",
        required=False,
        type=str,
        dest="temp_path",
        default="",
        help="Temp folder for chunk creation specification. Useful when using a cluster with a scratch folder",
    )

    # Global group
    parser.add_argument(
        "-n",
        "--first_n",
        required=False,
        type=int,
        dest="first_n",
        default=float("inf"),
        help="Select N reads to run on instead of all.",
    )
    parser.add_argument(
        "-o",
        "--output",
        required=False,
        type=str,
        default="Results",
        dest="outfolder",
        help="Results will be written to this folder",
    )
    parser.add_argument(
        "--dense",
        required=False,
        action="store_true",
        default=False,
        dest="dense",
        help="Add a dense output to the results folder",
    )
    parser.add_argument(
        "-u",
        "--unmapped-tags",
        required=False,
        type=str,
        dest="unmapped_file",
        default="unmapped.csv",
        help="Write table of unknown TAGs to file.",
    )
    parser.add_argument(
        "-ut",
        "--unknown-top-tags",
        required=False,
        dest="unknowns_top",
        type=int,
        default=100,
        help="Top n unmapped TAGs.",
    )
    parser.add_argument(
        "--debug", action="store_true", help="Print extra information for debugging."
    )
    parser.add_argument(
        "--version",
        action="version",
        version="CITE-seq-Count v{}".format(get_package_version()),
        help="Print version number.",
    )
    # Finally! Too many options XD
    return parser
