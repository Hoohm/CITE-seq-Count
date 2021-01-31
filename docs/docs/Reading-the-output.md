Mtx format
------------------------

The mtx, [matrix market](http://networkrepository.com/mtx-matrix-market-format.html), format is a sparse format for matrices. It only stores non zero values and is becoming popular in single-cell softwares.

The main advantage is that it requires less space than a dense matrix and that you can easily add different feature names within the same object.

For CITE-seq-Count, the output looks like this:

```
OUTFOLDER/
-- umi_count/
-- -- matrix.mtx.gz
-- -- features.tsv.gz
-- -- barcodes.tsv.gz
-- read_count/
-- -- matrix.mtx.gz
-- -- features.tsv.gz
-- -- barcodes.tsv.gz
-- unmapped.csv
-- run_report.yaml
```

File descriptions
-------------------

* `features.tsv.gz` contains the feature names, in this context our tags.
* `barcodes.tsv.gz` contains the cell barcodes. If running with a translation, first column is the translated barcode, second column is the original barcode found in the data.
* `matrix.mtx.gz` contains the actual values.
read_count and umi_count contain respectively the read counts and the collapsed umi counts. For analysis you should use the umi data. The read_count can be used to check if you have an overamplification or oversequencing issue with your protocol.
* `unmapped.csv` contains the top N tags that haven't been mapped.
* `run_report.yaml` contains the parameters used for the run as well as some statistics.
here is an example:

```
Date: 2019-10-01
Running time: 13.86 seconds
CITE-seq-Count Version: 1.4.3
Reads processed: 1000000
Percentage mapped: 33
Percentage unmapped: 67
Percentage too short: 0
  R1_too_short: 0
  R2_too_short: 0
Uncorrected cells: 0
Correction:
	Cell barcodes collapsing threshold: 1
	Cell barcodes corrected: 57
	UMI collapsing threshold: 2
	UMIs corrected: 329
Run parameters:
	Read1_filename: fastq/test_R1.fastq.gz,fastq/test2_R1.fastq.gz
	Read2_filename: fastq/test_R2.fastq.gz,fastq/test2_R2.fastq.gz
	Cell barcode:
		First position: 1
		Last position: 16
	UMI barcode:
		First position: 17
		Last position: 26
	Expected cells: 100
	Tags max errors: 1
	Start trim: 0
```

Packages to read MTX
--------------------------
**R**

I recommend using `Seurat` and their `Read10x` function to read the results.


With Seurat V3:

`Read10x('OUTFOLDER/umi_count/', gene.column=1)`

Version 1.5.0 of CSC came with some breaking changes. Older versions would use this command:

`Read10x('OUTFOLDER/umi_count/', gene.column=1)`

With Matrix:

```
library(Matrix)
matrix_dir = "/path_to_your_directory/out_cite_seq_count/umi_count/"
barcode.path <- paste0(matrix_dir, "barcodes.tsv.gz")
features.path <- paste0(matrix_dir, "features.tsv.gz")
matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")
mat <- readMM(file = matrix.path)
feature.names = read.delim(features.path, header = FALSE, stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path, header = FALSE, stringsAsFactors = FALSE)
colnames(mat) = barcode.names$V1
rownames(mat) = feature.names$V2
# rownames(mat) = feature.names$V1 if you are using an older version than 1.5.0
```

**Python**

I recommend using `scanpy` and their read_mtx function to read the results.

Example:

```
import scanpy
import pandas as pd
import os
path = 'umi_cell_corrected'
data = scanpy.read_mtx(os.path.join(path,'umi_count/matrix.mtx.gz'))
data = data.T
features = pd.read_csv(os.path.join(path, 'umi_count/features.tsv.gz'), header=None)
barcodes = pd.read_csv(os.path.join(path, 'umi_count/barcodes.tsv.gz'), header=None)
data.var_names = features[0]
data.obs_names = barcodes[1]
#data.obs_names = barcodes[0] if you are using an older version than 1.5.0
```
