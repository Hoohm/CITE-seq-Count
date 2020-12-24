import pytest
import random
import copy
from collections import Counter, namedtuple
from cite_seq_count import processing
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
    for i in range(n):
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
    import json
    from collections import defaultdict
    from collections import OrderedDict
    from collections import Counter

    from itertools import islice

    # Test file paths
    pytest.correct_R1_path = "tests/test_data/fastq/correct_R1.fastq.gz"
    pytest.correct_R2_path = "tests/test_data/fastq/correct_R2.fastq.gz"
    pytest.file_path = "tests/test_data/fastq/test_csv.csv"

    pytest.chunk_size = 800
    tag = namedtuple("tag", ["name", "sequence", "id"])
    pytest.tags = [
        tag(name="test1", sequence="CGTACGTAGCCTAGC", id=0),
        tag(name="test2", sequence="CGTAGCTCG", id=1),
        tag(name="unmapped", sequence="UNKNOWN", id=3),
    ]

    pytest.barcode_slice = slice(0, 16)
    pytest.umi_slice = slice(16, 26)
    pytest.correct_whitelist = set(["ACTGTTTTATTGGCCT", "TTCATAAGGTAGGGAT"])
    pytest.legacy = False
    pytest.debug = False
    pytest.start_trim = 0
    pytest.maximum_distance = 5
    pytest.results = {
        "ACTGTTTTATTGGCCT": {
            0: Counter({b"CATTAGTGGT": 3, b"CATTAGTGGG": 2, b"CATTCGTGGT": 1})
        },
        "TTCATAAGGTAGGGAT": {
            1: Counter({b"TAGCTTAGTA": 3, b"TAGCTTAGTC": 2, b"GCGATGCATA": 1})
        },
    }
    pytest.corrected_results = {
        "ACTGTTTTATTGGCCT": {0: Counter({b"CATTAGTGGT": 6})},
        "TTCATAAGGTAGGGAT": {1: Counter({b"TAGCTTAGTA": 5, b"GCGATGCATA": 1})},
    }
    pytest.umis_per_cell = Counter({"ACTGTTTTATTGGCCT": 1, "TTCATAAGGTAGGGAT": 2})
    pytest.reads_per_cell = Counter({"ACTGTTTTATTGGCCT": 3, "TTCATAAGGTAGGGAT": 6})
    pytest.expected_cells = 2
    pytest.no_match = Counter()
    pytest.collapsing_threshold = 1
    pytest.sliding_window = False
    pytest.max_umis = 20000

    pytest.sequence_pool = []
    pytest.tags_tuple = preprocessing.check_tags(
        preprocessing.parse_tags_csv("tests/test_data/tags/pass/correct.csv"), 5
    )[0]
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
            counts[processing.find_best_match(seq, pytest.tags_tuple, distance)] += 1
        assert counts[tag.id] == 4


@pytest.mark.dependency()
def test_find_best_match_with_2_distance(data):
    distance = 2
    for tag in pytest.tags_tuple:
        counts = Counter()
        if tag.name == "unmapped":
            continue
        for seq in extend_seq_pool(tag.sequence, distance):
            counts[processing.find_best_match(seq, pytest.tags_tuple, distance)] += 1
        assert counts[tag.id] == 4


@pytest.mark.dependency()
def test_find_best_match_with_3_distance(data):
    distance = 3
    for tag in pytest.tags_tuple:
        counts = Counter()
        if tag.name == "unmapped":
            continue
        for seq in extend_seq_pool(tag.sequence, distance):
            counts[processing.find_best_match(seq, pytest.tags_tuple, distance)] += 1
        assert counts[tag.id] == 4


@pytest.mark.dependency()
def test_find_best_match_with_3_distance_reverse(data):
    distance = 3
    for tag in pytest.tags_tuple:
        counts = Counter()
        if tag.name == "unmapped":
            continue
        for seq in extend_seq_pool(tag.sequence, distance):
            counts[processing.find_best_match(seq, pytest.tags_tuple, distance)] += 1
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
    (results, _) = processing.map_reads(pytest.mappint_input_test)
    print(results)
    assert len(results) == 2


@pytest.mark.dependency(depends=["test_classify_reads_multi_process"])
def test_correct_umis(data):
    temp = processing.correct_umis((pytest.results, 2, pytest.max_umis))
    results = temp[0]
    n_corrected = temp[1]
    for cell_barcode in results.keys():
        for TAG in results[cell_barcode]:
            assert len(results[cell_barcode][TAG]) == len(
                pytest.corrected_results[cell_barcode][TAG]
            )
            assert sum(results[cell_barcode][TAG].values()) == sum(
                pytest.corrected_results[cell_barcode][TAG].values()
            )
    assert n_corrected == 3


@pytest.mark.dependency(depends=["test_correct_umis"])
def test_correct_cells(data):
    processing.correct_cells(
        pytest.corrected_results,
        pytest.reads_per_cell,
        pytest.umis_per_cell,
        pytest.expected_cells,
        pytest.collapsing_threshold,
        pytest.tags,
    )


@pytest.mark.dependency(depends=["test_correct_umis"])
def test_generate_sparse_umi_matrices(data):
    umi_results_matrix = processing.generate_sparse_matrices(
        pytest.corrected_results,
        pytest.tags,
        set(["ACTGTTTTATTGGCCT", "TTCATAAGGTAGGGAT"]),
        umi_counts=True,
    )
    assert umi_results_matrix.shape == (3, 2)
    total_umis = 0
    for i in range(umi_results_matrix.shape[0]):
        for j in range(umi_results_matrix.shape[1]):
            total_umis += umi_results_matrix[i, j]
    assert total_umis == 3


@pytest.mark.dependency(depends=["test_correct_umis"])
def test_generate_sparse_read_matrices(data):
    read_results_matrix = processing.generate_sparse_matrices(
        pytest.corrected_results,
        pytest.tags,
        set(["ACTGTTTTATTGGCCT", "TTCATAAGGTAGGGAT"]),
        umi_counts=False,
    )
    assert read_results_matrix.shape == (3, 2)
    total_umis = 0
    for i in range(read_results_matrix.shape[0]):
        for j in range(read_results_matrix.shape[1]):
            total_umis += read_results_matrix[i, j]
    assert total_umis == 12
