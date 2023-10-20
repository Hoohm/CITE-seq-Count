import pytest
import random
import copy
from collections import Counter, namedtuple
from cite_seq_count import mapping
from cite_seq_count import preprocessing


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
def data():
    from collections import Counter

    pytest.file_path = "tests/test_data/fastq/test_csv.csv"
    pytest.debug = False
    pytest.barcode_slice = slice(0, 16)
    pytest.umi_slice = slice(16, 26)
    pytest.correct_reference_list = set(["ACTGTTTTATTGGCCT", "TTCATAAGGTAGGGAT"])
    pytest.maximum_distance = 5
    pytest.results = {
        "ACTGTTTTATTGGCCT": {
            0: Counter({b"CATTAGTGGT": 3, b"CATTAGTGGG": 2, b"CATTCGTGGT": 1})
        },
        "TTCATAAGGTAGGGAT": {
            1: Counter({b"TAGCTTAGTA": 3, b"TAGCTTAGTC": 2, b"GCGATGCATA": 1})
        },
    }

    pytest.sliding_window = False
    pytest.sequence_pool = []
    pytest.tags_tuple = preprocessing.parse_tags_csv(
        preprocessing.parse_tags_csv("tests/test_data/tags/pass/correct.csv")
    )
    pytest.mapping_input = namedtuple(
        "mapping_input",
        ["filename", "tags", "debug", "maximum_distance", "sliding_window"],
    )
    pytest.mappint_input_test = pytest.mapping_input(
        filename=pytest.file_path,
        tags=pytest.tags_tuple,
        debug=pytest.debug,
        maximum_distance=pytest.maximum_distance,
        sliding_window=pytest.sliding_window,
    )


@pytest.mark.dependency()
def test_find_best_match_with_1_distance(data):
    distance = 1
    for tag in pytest.tags_tuple:
        counts = Counter()
        if tag.name == "unmapped":
            continue
        for seq in extend_seq_pool(tag.sequence, distance):
            counts[mapping.find_best_match(seq, pytest.tags_tuple, distance)] += 1
        assert counts[tag.id] == 4


@pytest.mark.dependency()
def test_find_best_match_with_2_distance(data):
    distance = 2
    for tag in pytest.tags_tuple:
        counts = Counter()
        if tag.name == "unmapped":
            continue
        for seq in extend_seq_pool(tag.sequence, distance):
            counts[mapping.find_best_match(seq, pytest.tags_tuple, distance)] += 1
        assert counts[tag.id] == 4


@pytest.mark.dependency()
def test_find_best_match_with_3_distance(data):
    distance = 3
    for tag in pytest.tags_tuple:
        counts = Counter()
        for seq in extend_seq_pool(tag.sequence, distance):
            counts[mapping.find_best_match(seq, pytest.tags_tuple, distance)] += 1
        assert counts[tag.id] == 4


@pytest.mark.dependency()
def test_find_best_match_with_3_distance_reverse(data):
    distance = 3
    for tag in pytest.tags_tuple:
        counts = Counter()
        if tag.name == "unmapped":
            continue
        for seq in extend_seq_pool(tag.sequence, distance):
            counts[mapping.find_best_match(seq, pytest.tags_tuple, distance)] += 1
        assert counts[tag.id] == 4


@pytest.mark.dependency(
    depends=[
        "test_find_best_match_with_1_distance",
        "test_find_best_match_with_2_distance",
        "test_find_best_match_with_3_distance",
        "test_find_best_match_with_3_distance_reverse",
    ]
)
def test_classify_reads_multi_process(data):
    (results, _) = mapping.map_reads(pytest.mappint_input_test)
    print(results)
    assert len(results) == 2
