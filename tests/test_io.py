import pytest
import os
import gzip
from pathlib import Path
from polars.testing import assert_frame_equal
from cite_seq_count.constants import (
    BARCODE_COLUMN,
    COUNT_COLUMN,
    FEATURE_NAME_COLUMN,
    SEQUENCE_COLUMN,
    SUBSET_COLUMN,
    FEATURE_ID_COLUMN,
    BARCODE_ID_COLUMN,
    UNMAPPED_NAME,
)
from cite_seq_count.io import (
    get_n_lines,
    get_read_paths,
    write_data_to_mtx,
    create_mtx_df,
    write_mapping_input_from_fastqs,
    read_R1_polars,
    read_R2_polars,
    write_fastq_inputs_as_parquet,
    concat_reads,
)

from cite_seq_count.chemistry import Chemistry
import polars as pl

# copied from https://stackoverflow.com/questions/3431825/generating-an-md5-checksum-of-a-file
import hashlib


def md5(fname):
    hash_md5 = hashlib.md5()
    if fname.endswith("gz"):
        with gzip.open(fname, "rb") as f:
            string = f.read()
            hash_md5.update(string)
    else:
        with open(fname, "r") as f:
            string = f.read()
            hash_md5.update(string.encode())
    return hash_md5.hexdigest()


@pytest.fixture
def chemistry_def():
    return Chemistry(
        name="test",
        cell_barcode_start=1,
        cell_barcode_end=16,
        umi_barcode_start=17,
        umi_barcode_end=26,
        r2_trim_start=0,
        barcode_reference_path="tests/reference_lists/pass/translation.csv",
    )


@pytest.fixture
def correct_R1():
    return Path("tests/test_data/fastq/correct_R1.fastq.gz")


@pytest.fixture
def correct_R2():
    return Path("tests/test_data/fastq/correct_R2.fastq.gz")


@pytest.fixture
def corrupt_R1():
    return Path("tests/test_data/fastq/corrupted_R1.fastq.gz")


@pytest.fixture
def corrupt_R2():
    return Path("tests/test_data/fastq/corrupted_R2.fastq.gz")


@pytest.fixture
def correct_R1_multi():
    return "tests/test_data/fastq/correct_R1.fastq.gz,tests/test_data/fastq/correct_R1.fastq.gz"


@pytest.fixture
def correct_R2_multi():
    return "tests/test_data/fastq/correct_R2.fastq.gz,tests/test_data/fastq/correct_R2.fastq.gz"


@pytest.fixture
def corrupt_R2_multi():
    return "tests/test_data/fastq/correct_R2.fastq.gz,tests/test_data/fastq/corrupted_R2.fastq.gz"


@pytest.fixture
def incorrect_R2_multipath():
    return "path/to/R2_1.fastq.gz,path/to/R2_2.fastq.gz,path/to/R2_3.fastq.gz"


@pytest.fixture
def correct_multipath_combined(correct_R1, correct_R2):
    return (
        [correct_R1, correct_R1],
        [correct_R2, correct_R2],
    )


@pytest.mark.dependency()
def test_get_n_lines(correct_R1):
    assert get_n_lines(correct_R1) == (200 * 4)


def test_corrrect_multipath(
    correct_R1_multi, correct_R2_multi, correct_multipath_combined
):
    assert (
        get_read_paths(correct_R1_multi, correct_R2_multi) == correct_multipath_combined
    )


@pytest.mark.dependency(depends=["test_get_n_lines"])
def test_incorrect_multipath(correct_R1_multi, incorrect_R2_multipath):
    with pytest.raises(SystemExit):
        get_read_paths(correct_R1_multi, incorrect_R2_multipath)


@pytest.mark.dependency(depends=["test_get_n_lines"])
def test_get_n_lines_not_multiple_of_4(corrupt_R1):
    with pytest.raises(SystemExit):
        get_n_lines(corrupt_R1)


def test_write_data_to_mtx(tmp_path):
    # Create test data
    main_df = pl.LazyFrame(
        {
            FEATURE_NAME_COLUMN: ["test1", "test2", "test3"],
            BARCODE_COLUMN: [
                "ACTGTTTTATTGGCCT",
                "TTCATAAGGTAGGGAT",
                "ACTGTTTTATTGGCCT",
            ],
            COUNT_COLUMN: [10, 20, 30],
        }
    )
    tags_df = pl.DataFrame(
        {
            FEATURE_NAME_COLUMN: ["test1", "test2", "test3"],
            SEQUENCE_COLUMN: ["CGTA", "CGTA", "CGTA"],
        }
    )
    subset_df = pl.LazyFrame({SUBSET_COLUMN: ["ACTGTTTTATTGGCCT", "TTCATAAGGTAGGGAT"]})
    data_type = "umi"
    outpath = str(tmp_path)
    print(outpath)

    # Call the function
    write_data_to_mtx(
        main_df=main_df,
        tags_df=tags_df,
        subset_df=subset_df,
        data_type=data_type,
        outpath=outpath,
    )

    # Check if the output files exist
    assert os.path.exists(os.path.join(outpath, "umi_count", "matrix.mtx.gz"))
    assert os.path.exists(os.path.join(outpath, "umi_count", "features.tsv.gz"))
    assert os.path.exists(os.path.join(outpath, "umi_count", "barcodes.tsv.gz"))

    # TODO: Add assertions to validate the content of the output files


def test_create_mtx_df():
    # Create test data
    main_df = pl.LazyFrame(
        {
            FEATURE_NAME_COLUMN: ["test1", "test2", "test3"],
            BARCODE_COLUMN: [
                "ACTGTTTTATTGGCCT",
                "TTCATAAGGTAGGGAT",
                "ACTGTTTTATTGGCCT",
            ],
            COUNT_COLUMN: [10, 20, 30],
        }
    )
    tags_df = pl.DataFrame(
        {
            FEATURE_NAME_COLUMN: ["test1", "test2", "test3"],
            SEQUENCE_COLUMN: ["CGTA", "CGTA", "CGTA"],
        }
    )
    subset_df = pl.LazyFrame({SUBSET_COLUMN: ["ACTGTTTTATTGGCCT", "TTCATAAGGTAGGGAT"]})

    # Call the function
    mtx_df, tags_indexed, barcodes_indexed = create_mtx_df(main_df, tags_df, subset_df)

    # Check the output
    expected_mtx_df = pl.DataFrame(
        {
            FEATURE_ID_COLUMN: [1, 2, 3],
            BARCODE_ID_COLUMN: [1, 2, 1],
            COUNT_COLUMN: [10, 20, 30],
        },
        schema={
            FEATURE_ID_COLUMN: pl.UInt32,
            BARCODE_ID_COLUMN: pl.UInt32,
            COUNT_COLUMN: pl.Int64,
        },
    )
    expected_tags_indexed = pl.DataFrame(
        {
            FEATURE_NAME_COLUMN: ["test1", "test2", "test3", UNMAPPED_NAME],
            SEQUENCE_COLUMN: ["CGTA", "CGTA", "CGTA", "UNKNOWN"],
            FEATURE_ID_COLUMN: [1, 2, 3, 4],
        },
        schema={
            FEATURE_ID_COLUMN: pl.UInt32,
            SEQUENCE_COLUMN: pl.String,
            FEATURE_NAME_COLUMN: pl.String,
        },
    )
    expected_barcodes_indexed = pl.DataFrame(
        {
            SUBSET_COLUMN: ["ACTGTTTTATTGGCCT", "TTCATAAGGTAGGGAT"],
            BARCODE_ID_COLUMN: [1, 2],
        },
        schema={SUBSET_COLUMN: pl.String, BARCODE_ID_COLUMN: pl.UInt32},
    )

    assert_frame_equal(mtx_df, expected_mtx_df)
    assert_frame_equal(tags_indexed, expected_tags_indexed, check_column_order=False)
    assert_frame_equal(
        barcodes_indexed, expected_barcodes_indexed, check_column_order=False
    )


@pytest.fixture
def fastq_reads():
    return pl.scan_csv(
        "tests/test_data/fastq/test_csv.csv",
        has_header=False,
        schema={"barcode": pl.String, "umi": pl.String, "r2": pl.String},
    )


@pytest.fixture
def fastq_paths(correct_R1, correct_R2):
    return [
        (correct_R1, correct_R2),
    ]


@pytest.fixture
def correct_R1_csv():
    return pl.scan_csv("tests/test_data/fastq/correct_R1.csv", has_header=True)


@pytest.fixture
def correct_R2_csv():
    return pl.scan_csv("tests/test_data/fastq/correct_R2.csv", has_header=True)


def test_read_R1_polars(correct_R1, chemistry_def, correct_R1_csv):
    # Call the function
    result = read_R1_polars(correct_R1, chemistry_def)
    # Define the expected output
    # Check the output
    assert_frame_equal(result, correct_R1_csv)


def test_read_R2_polars(correct_R2, correct_R2_csv, chemistry_def):
    # Call the function
    result = read_R2_polars(correct_R2, r2_min_length=20, chemistry=chemistry_def)
    # Define the expected output

    # Check the output
    assert_frame_equal(result, correct_R2_csv)


def test_concat_reads(correct_R1, correct_R2, chemistry_def, tmp_path):
    # Call the function to be tested
    temp_file_path, r1_too_short, r2_too_short, total_reads = concat_reads(
        correct_R1, correct_R2, chemistry_def, 100, 1, 20, tmp_path
    )
    os.remove(temp_file_path)
    # Assert the expected output
    assert r1_too_short == 0
    assert r2_too_short == 0
    assert total_reads == 100


def test_write_fastq_inputs_as_parquet(tmp_path, correct_R1, correct_R2, chemistry_def):
    # Set up test data
    read_paths = [(correct_R1, correct_R2)]
    temp_path = tmp_path
    r2_min_length = 20
    top_n_reads = 100

    # Call the function
    r1_too_short, r2_too_short, total_reads, temp_file_list = write_fastq_inputs_as_parquet(
        read_paths=read_paths,
        temp_path=temp_path,
        chemistry_def=chemistry_def,
        r2_min_length=r2_min_length,
        top_n_reads=top_n_reads,
    )

    # Assert the result
    print(temp_file_list)
    assert r1_too_short == 0
    assert r2_too_short == 0
    assert total_reads == top_n_reads
    assert temp_path.exists()
