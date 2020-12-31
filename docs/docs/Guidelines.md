Guidelines for typical chemistries
------------------------

### 10x Genomics

## 10xV3

10x genomics V3 chemistry for feature barcoding is using a mapping between RNA cell barcodes and Protein cell barcodes.
You can find this mapping [here](https://github.com/10XGenomics/cellranger/blob/master/lib/python/cellranger/barcodes/translation/3M-february-2018.txt.gz)

## POST 1.5.0 instructions

Since version 1.5.0, this is taken care of by CSC if you provide the translation column in the reference `--reference_file` file as described in the documentation.
* The dense output will have the translated barcodes in the header.
* the MTX output will have two columns. The first column is the translated barcode given by the reference list csv and the second column is the original barcode found in the data.
* I recommend using the MTX format because it contains both cell barcodes.

## PRE 1.5.0 instructions

Since the list is composed of ~7M cells instead of the ~3M described in the technologie, using the reference_list of V3 as an input is unwise.
I suggest running CITE-seq-Count with using the `-cells` argument instead.