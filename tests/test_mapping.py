from ast import parse
import pytest
import random
import copy
from cite_seq_count import mapping
from polars.testing import assert_frame_equal
from cite_seq_count.preprocessing import parse_tags_csv
from cite_seq_count.constants import R2_COLUMN, SEQUENCE_COLUMN, FEATURE_NAME_COLUMN
import polars as pl


def complete_poly_A(seq, final_length=40):
    poly_A_len = final_length - len(seq)
    return seq + "A" * poly_A_len


def get_sequences(ref_path):
    sequences = []
    with open(ref_path, "r") as adt_ref:
        lines = adt_ref.readlines()
        entries = int(len(lines) / 2)
        for i in range(0, entries, 2):
            sequences.append(complete_poly_A(lines[i + 1].strip()))
    return sequences


def extend_seq_pool(ref_seq, distance):
    extended_pool = [complete_poly_A(ref_seq)]
    extended_pool.append(modify(ref_seq, distance, modification_type="mutate"))
    extended_pool.append(modify(ref_seq, distance, modification_type="mutate"))
    extended_pool.append(modify(ref_seq, distance, modification_type="mutate"))
    return extended_pool


def modify(seq, n, modification_type):
    bases = list("ATGCN")
    positions = list(range(len(seq)))
    seq = list(seq)
    for _ in range(n):
        if modification_type == "mutate":
            position = random.choice(positions)
            positions.remove(position)
            temp_bases = copy.copy(bases)
            del temp_bases[bases.index(seq[position])]
            seq[position] = random.choice(temp_bases)
        elif modification_type == "delete":
            del seq[random.randint(0, len(seq) - 2)]
        elif modification_type == "add":
            position = random.randint(0, len(seq) - 1)
            seq.insert(position, random.choice(bases))
    return complete_poly_A("".join(seq))


@pytest.fixture
def small_dataset_path():
    return "tests/test_data/fastq/test_csv.csv"


@pytest.fixture
def parsed_tags_df():
    return parse_tags_csv("tests/test_data/tags/pass/correct.csv")


@pytest.fixture
def r2_df():
    # Create a sample DataFrame for r2_df
    return pl.LazyFrame(
        {
            R2_COLUMN: ["ACTGTTTTATTGGCCT", "TTCATAAGGTAGGGAT", "AGCTAGCTAGCTAGCT"],
        }
    )


@pytest.fixture
def parsed_tags():
    # Create a sample DataFrame for parsed_tags
    return pl.DataFrame(
        {
            SEQUENCE_COLUMN: ["ACTGTTTTATTGGCCT", "TTCATAAGGTAGGGAT"],
            FEATURE_NAME_COLUMN: ["feature1", "feature2"],
        }
    )


def test_map_reads_polars_with_dist1(r2_df, parsed_tags):
    maximum_distance = 1
    expected_mapped = pl.LazyFrame(
        {
            R2_COLUMN: ["ACTGTTTTATTGGCCT", "TTCATAAGGTAGGGAT"],
            FEATURE_NAME_COLUMN: ["feature1", "feature2"],
        }
    )
    expected_unmapped = pl.LazyFrame(
        {R2_COLUMN: ["AGCTAGCTAGCTAGCT"], FEATURE_NAME_COLUMN: ["unmapped"]}
    )

    mapped, unmapped = mapping.map_reads_polars(r2_df, parsed_tags, maximum_distance)

    assert_frame_equal(mapped, expected_mapped)
    assert_frame_equal(unmapped, expected_unmapped)


def test_map_reads_polars_with_dist2(r2_df, parsed_tags):
    maximum_distance = 2
    expected_mapped = pl.LazyFrame(
        {
            R2_COLUMN: ["ACTGTTTTATTGGCCT", "TTCATAAGGTAGGGAT"],
            FEATURE_NAME_COLUMN: ["feature1", "feature2"],
        }
    )
    expected_unmapped = pl.LazyFrame(
        {R2_COLUMN: ["AGCTAGCTAGCTAGCT"], FEATURE_NAME_COLUMN: ["unmapped"]}
    )

    mapped, unmapped = mapping.map_reads_polars(r2_df, parsed_tags, maximum_distance)

    assert_frame_equal(mapped, expected_mapped)
    assert_frame_equal(unmapped, expected_unmapped)
