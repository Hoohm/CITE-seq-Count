# REQUIRED_TAGS_HEADER = ["sequence", "feature_name"]
REQUIRED_CELLS_REF_HEADER = ["reference"]
OPTIONAL_CELLS_REF_HEADER = ["translation"]
# Polars column names
# Tags input
FEATURE_NAME_COLUMN = "feature_name"
SEQUENCE_COLUMN = "sequence"
REQUIRED_TAGS_HEADER = [FEATURE_NAME_COLUMN, SEQUENCE_COLUMN]
UNMAPPED_NAME = "unmapped"

# Reads input
BARCODE_COLUMN = "barcode"
CORRECTED_BARCODE_COLUMN = "corrected_barcode"
UMI_COLUMN = "umi"
R2_COLUMN = "r2"
COUNT_COLUMN = "count"
# Barcode input
REFERENCE_COLUMN = "reference"
TRANSLATION_COLUMN = "translation"
SUBSET_COLUMN = "subset"
STRIP_CHARS = '"0123456789- \t\n'

# MTX format
BARCODE_ID_COLUMN = "barcode_id"
FEATURE_ID_COLUMN = "feature_id"
MTX_HEADER = """%%MatrixMarket matrix coordinate integer general\n%\n"""
TEMP_MTX = "temp.mtx"
FEATURES_MTX = "features.tsv.gz"
BARCODE_MTX = "barcodes.tsv.gz"
MATRIX_MTX = "matrix.mtx.gz"
