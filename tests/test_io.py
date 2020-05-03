import pytest
from cite_seq_count import io



@pytest.fixture
def data():
    from collections import OrderedDict
    from scipy import sparse
    test_matrix = sparse.dok_matrix((4,2))
    test_matrix[1,1] = 1
    pytest.sparse_matrix = test_matrix
    pytest.top_cells = set(['ACTGTTTTATTGGCCT','TTCATAAGGTAGGGAT'])
    pytest.ordered_tags_map = OrderedDict({
        'test3':{'id':0, 'sequence': 'CGTA'},
        'test2':{'id':1, 'sequence': 'CGTA'},
        'test1': {'id':3, 'sequence': 'CGTA'},
        'unmapped': {'id':4, 'sequence': 'CGTA'}
        })
    pytest.data_type = 'umi'
    pytest.outfolder = 'tests/test_data/'

def test_write_to_files(data, tmpdir):
    import gzip
    import scipy
    io.write_to_files(pytest.sparse_matrix,
        pytest.top_cells,
        pytest.ordered_tags_map,
        pytest.data_type,
        tmpdir)
    file = tmpdir.join('umi_count/matrix.mtx.gz')
    with gzip.open(file, 'rb') as mtx_file:
        assert isinstance(scipy.io.mmread(mtx_file) ,scipy.sparse.coo.coo_matrix)
