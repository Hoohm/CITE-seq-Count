# Change Log
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/)
and this project adheres to [Semantic Versioning](http://semver.org/).

## [1.3.1a] - 2018-10-22
### Added
 - Printing version now from the help
 - Added a possibility to store unmatched tags in a file using the option `-u`


## [1.3] - 2018-09-05
### Added
- `-l` `--legacy` option now available. This option will then end of the regex to find the TAG barcode in R2. Use `-l` if you have TAG barcodes ending with [TGC] + polyA tails.
- New warning that checks R1 length in the fastq file and the cell and umi barcodes given as input.

### Changed
- TAG barcode sequence are now added to the name of the tag in the rows of the results. This will help reproducibility since the barcode sequence will already be in the count results.


## [1.2] - 2018-08-02
### Added
- pandas dependcy


## [1.1] - 2018-08-02
### Changed
- Cite-seq-Count is now a python package!

## [0.2] - 2018-07-17
### Added
- Compatibility with tags of different lengths
- Possibility to use a whitelist of barcodes to extract
- Uses now fuzzy regex for structure detection

### Changed
- Regex is now optional.
- You can now use only one tag, the script won't crash


## [0.1] - 2017-08-07
### Added
- Processing of CITE-seq data with antibody tags and/or HTO