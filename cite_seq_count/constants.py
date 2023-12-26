

# REQUIRED_TAGS_HEADER = ["sequence", "feature_name"]
REQUIRED_CELLS_REF_HEADER = ["reference"]
OPTIONAL_CELLS_REF_HEADER = ["translation"]
# Polars column names
# Tags input
FEATURE_NAME_COLUMN = "feature_name"
SEQUENCE_COLUMN = "sequence"
REQUIRED_TAGS_HEADER = [FEATURE_NAME_COLUMN, SEQUENCE_COLUMN]
# Reads input
BARCODE_COLUMN = "barcode"
CORRECTED_BARCODE_COLUMN = "corrected_barcode"
UMI_COLUMN = "umi"
R2_COLUMN = "r2"
# Barcode input
REFERENCE_COLUMN = "reference"
TRANSLATION_COLUMN = "translation"
WHITELIST_COLUMN = "whitelist"
STRIP_CHARS = '"0123456789- \t\n'

UNMAPPED_NAME = "unmapped"
