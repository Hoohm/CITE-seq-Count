import pytest
from collections import namedtuple
from cite_seq_count import processing


@pytest.fixture
def data():
    from collections import Counter

    tag = namedtuple("tag", ["name", "sequence", "id"])
    pytest.tags = [
        tag(name="test1", sequence="CGTACGTAGCCTAGC", id=0),
        tag(name="test2", sequence="CGTAGCTCG", id=1),
    ]
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
    pytest.collapsing_threshold = 1
    pytest.max_umis = 20000


@pytest.mark.dependency()
def test_correct_umis(data):
    temp = processing.correct_umis_in_cells((pytest.results, 2, pytest.max_umis, 2))
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
    processing.correct_cells_no_translation_list(
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
        ["ACTGTTTTATTGGCCT", "TTCATAAGGTAGGGAT"],
        umi_counts=True,
    )
    assert umi_results_matrix.shape == (2, 2)
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
        ["ACTGTTTTATTGGCCT", "TTCATAAGGTAGGGAT"],
        umi_counts=False,
    )
    assert read_results_matrix.shape == (3, 2)
    total_umis = 0
    for i in range(read_results_matrix.shape[0]):
        for j in range(read_results_matrix.shape[1]):
            total_umis += read_results_matrix[i, j]
    assert total_umis == 12
