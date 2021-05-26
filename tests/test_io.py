import pytest
import os
import gzip
import scipy
from cite_seq_count import io
from collections import namedtuple
import numpy as np

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
def data():
    from collections import OrderedDict
    from scipy import sparse

    test_matrix = sparse.dok_matrix((4, 2), dtype=np.int32)
    test_matrix[1, 1] = 1
    pytest.sparse_matrix = test_matrix
    pytest.filtered_cells = ["ACTGTTTTATTGGCCT", "TTCATAAGGTAGGGAT"]
    tag = namedtuple("tag", ["name", "sequence", "id"])
    pytest.ordered_tags_map = [
        tag(name="test1", sequence="CGTA", id=0),
        tag(name="test2", sequence="CGTA", id=1),
        tag(name="test3", sequence="CGTA", id=2),
        tag(name="unmapped", sequence="UNKNOWN", id=3),
    ]

    pytest.data_type = "umi"


def test_write_to_files_wo_translation(data, tmpdir):
    reference_dict = {"ACTGTTTTATTGGCCT": 0, "TTCATAAGGTAGGGAT": 0}
    output_path = os.path.join(tmpdir, "without_translation")

    mtx_path = os.path.join(output_path, "umi_count", "matrix.mtx.gz")
    features_path = os.path.join(output_path, "umi_count", "features.tsv.gz")
    barcodes_path = os.path.join(output_path, "umi_count", "barcodes.tsv.gz")
    md5_sums = {
        barcodes_path: "b7af6a32e83963606f181509a571966f",
        features_path: "e889e780dbce481287c993dd043714c8",
        mtx_path: "3ea98c44d88a947215bace0c72ac1303",
    }

    io.write_to_files(
        pytest.sparse_matrix,
        pytest.filtered_cells,
        pytest.ordered_tags_map,
        pytest.data_type,
        output_path,
        reference_dict=reference_dict,
    )
    file_path = os.path.join(tmpdir, "without_translation", "umi_count/matrix.mtx.gz")
    with gzip.open(file_path, "rb") as mtx_file:
        assert isinstance(scipy.io.mmread(mtx_file), scipy.sparse.coo.coo_matrix)
    assert md5_sums[barcodes_path] == md5(barcodes_path)
    assert md5_sums[features_path] == md5(features_path)
    assert md5_sums[mtx_path] == md5(mtx_path)


def test_write_to_files_with_translation(data, tmpdir):
    reference_dict = {
        "ACTGTTTTATTGGCCT": "GGCTTCGATACTAGAT",
        "TTCATAAGGTAGGGAT": "GATCGGATAGCTAATA",
    }
    output_path = os.path.join(tmpdir, "with_translation")

    mtx_path = os.path.join(output_path, "umi_count", "matrix.mtx.gz")
    features_path = os.path.join(output_path, "umi_count", "features.tsv.gz")
    barcodes_path = os.path.join(output_path, "umi_count", "barcodes.tsv.gz")

    md5_sums = {
        barcodes_path: "fce83378b4dd548882fb9271bdd5b4f1",
        features_path: "e889e780dbce481287c993dd043714c8",
        mtx_path: "3ea98c44d88a947215bace0c72ac1303",
    }

    io.write_to_files(
        pytest.sparse_matrix,
        pytest.filtered_cells,
        pytest.ordered_tags_map,
        pytest.data_type,
        output_path,
        reference_dict=reference_dict,
    )
    file_path = os.path.join(tmpdir, "with_translation", "umi_count/matrix.mtx.gz")
    with gzip.open(file_path, "rb") as mtx_file:
        assert isinstance(scipy.io.mmread(mtx_file), scipy.sparse.coo.coo_matrix)
    assert md5_sums[barcodes_path] == md5(barcodes_path)
    assert md5_sums[features_path] == md5(features_path)
    assert md5_sums[mtx_path] == md5(mtx_path)


def test_write_to_dense_wo_translation(data, tmpdir):
    reference_dict = {"ACTGTTTTATTGGCCT": 0, "TTCATAAGGTAGGGAT": 0}
    output_path = os.path.join(tmpdir, "without_translation")
    csv_name = "dense_umis.tsv"
    csv_path = os.path.join(output_path, csv_name)

    md5_sums = {
        csv_path: "fef502237900ec386d100169fa1fab7c",
    }

    io.write_dense(
        sparse_matrix=pytest.sparse_matrix,
        ordered_tags=pytest.ordered_tags_map,
        columns=pytest.filtered_cells,
        outfolder=output_path,
        filename=csv_name,
    )
    file_path = os.path.join(tmpdir, "without_translation", csv_name)
    assert md5_sums[csv_path] == md5(file_path)
