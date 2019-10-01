Mtx format
------------------------

The mtx, [matrix market](http://networkrepository.com/mtx-matrix-market-format.html), format is a sparse format for matrices. It only stores non zero values and is becoming popular in single-cell softwares.

The main advantage is that it requires less space than a dense matrix and that you can easily add different feature names withint the same object.

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
* `barcodes.tsv.gz` contains the cell barcodes.
* `matrix.mtx.gz` contains the actual values.
read_count and umi_count contain respectively the read counts and the collapsed umi counts. For analysis you should use the umi data. The read_count can be used to check if you have an overamplification or oversequencing issue with your protocol.
* `unmapped.csv` contains the top N tags that haven't been mapped.
* `run_report.yaml` contains the parameters used for the run as well as some statistics.
her is an example:

```
Date: 2019-01-01
Running time: 50 minutes, 30 seconds
CITE-seq-Count Version: 1.4.0
Reads processed: 50000000
Percentage mapped: 95
Percentage unmapped: 5
Parameters:
	Read1_filename: read1.fastq.gz
	Read2_filename: read1.fastq.gz
	Cell barcode:
		First position: 1
		Last position: 16
	UMI barcode:
		First position: 17
		Last position: 26
	Tags max errors: 3
Correction:
	Cell barcodes collapsing threshold: 1
	Cell barcodes corrected: 20000
	UMI collapsing threshold: 2
	UMIs corrected: 30000

```

Packages to read MTX
--------------------------
**R**

I recommend using [Seurat](https://satijalab.org/seurat/) and their `Read10x` function to read the results.

Example: `Read10x('OUTFOLDER/umi_count', gene.column=1)`

**Python**

I recommend using [Scanpy](https://icb-scanpy.readthedocs-hosted.com/en/latest/index.html) and their read_mtx function to read the results.

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
data.obs_names = barcodes[0]
```
