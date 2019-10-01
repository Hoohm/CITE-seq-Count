Guidelines for typical chemistries
------------------------

### 10x Genomics

## 10xV3

10x genomics V3 chemistry for feature barcoding is using a mapping between RNA cell barcodes and Protein cell barcodes.
You can find this mapping [here](https://github.com/10XGenomics/cellranger/blob/master/lib/python/cellranger/barcodes/translation/3M-february-2018.txt.gz)

Since the list is composed of ~7M cells instead of the ~3M described in the technologie, using the whitelist of V3 as an input is unwise.
I suggest running CITE-seq-Count with using the `-cells` argument instead.