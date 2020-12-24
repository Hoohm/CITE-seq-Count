import pytest
from cite_seq_count import io
from collections import namedtuple


@pytest.fixture
def data():
    from collections import OrderedDict
    from scipy import sparse

    test_matrix = sparse.dok_matrix((4, 2))
    test_matrix[1, 1] = 1
    pytest.sparse_matrix = test_matrix
    pytest.filtered_cells = set(["ACTGTTTTATTGGCCT", "TTCATAAGGTAGGGAT"])
    tag = namedtuple("tag", ["name", "sequence", "id"])
    pytest.ordered_tags_map = [
        tag(name="test1", sequence="CGTA", id=0),
        tag(name="test2", sequence="CGTA", id=1),
        tag(name="test3", sequence="CGTA", id=2),
        tag(name="unmapped", sequence="UNKNOWN", id=3),
    ]

    pytest.data_type = "umi"
    pytest.outfolder = "tests/test_data/"


def test_write_to_files(data, tmpdir):
    import gzip
    import scipy

    io.write_to_files(
        pytest.sparse_matrix,
        pytest.filtered_cells,
        pytest.ordered_tags_map,
        pytest.data_type,
        tmpdir,
    )
    file = tmpdir.join("umi_count/matrix.mtx.gz")
    with gzip.open(file, "rb") as mtx_file:
        assert isinstance(scipy.io.mmread(mtx_file), scipy.sparse.coo.coo_matrix)
