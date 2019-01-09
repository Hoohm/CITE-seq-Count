import os
import gzip
import shutil

from scipy import io

def write_to_files(sparse_matrix, top_cells, ordered_tags_map, data_type, outfile):
    """Write the umi and read sparse matrices to file in gzipped mtx format.

    Args:
        sparse_matrix (dok_matrix): Results in a sparse matrix.
        final_results (dict): Results in a dict of dicts of Counters.
        ordered_tags_map (dict): Tags in order with indexes as values.
        data_type (string): A string definning if the data is umi or read based.
        oufile (string): Path to the mtx file.
    
    """
    prefix = data_type + '_count'
    os.makedirs(prefix, exist_ok=True)
    io.mmwrite(os.path.join(prefix,'matrix.mtx'),sparse_matrix)
    with gzip.open(os.path.join(prefix,'barcodes.tsv.gz'), 'wb') as barcode_file:
        for barcode in top_cells:
            barcode_file.write('{}\n'.format(barcode).encode())
    with gzip.open(os.path.join(prefix,'features.tsv.gz'), 'wb') as feature_file:
        for feature in ordered_tags_map:
            feature_file.write('{}\n'.format(feature).encode())
    with open(os.path.join(prefix,'matrix.mtx'),'rb') as mtx_in:
        with gzip.open(os.path.join(prefix,'matrix.mtx') + '.gz','wb') as mtx_gz:
            shutil.copyfileobj(mtx_in, mtx_gz)
    os.remove(os.path.join(prefix,'matrix.mtx'))