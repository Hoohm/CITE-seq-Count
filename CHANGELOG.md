# Change Log
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/)
and this project adheres to [Semantic Versioning](http://semver.org/).

## [1.4.3] - XX.03.2019
### Added
  - Support for multiple files as input. This allows you to not merge different lanes before
    running CITE-seq-Count. Thanks to @arkal for the contribution on this! #79

### Changed
  - Cell barcode whitelist is now based on reads instead of UMIs. Fixing issue #35,#37,#48,#49,#69
  - Added a low filter for UMI correction not trying to correct unique UMIs. Small improvement on performance.
  - Corrected typo from `--no_umi_correctio` to `--no_umi_correction`. Thanks for the catch @ChristophH. #72

## [1.4.2] - 11.03.2019
### Added
- Cell barcode correction based on whitelist instead of top cells. This will trigger automatically
  when a whitelist is provided.
- UMI correction for tags with more than than 20'000 umis is too slow right now.
  Cells containing a TAG with more than 20'000 umis will be discarded. This means that more than 20'000
  antibody tags with different umis have been captured which is odd. Such high numbers could mean
  that cells have been aggregating and can't be used for downstream analysis. 
  The number of cells are reported under `Bad cells` and the dense matrix of those cells is provided in `uncorrected_cells`.
  Fixing issues #32
- The option to not correct for umis with `--no_umi_correction`.
- Now prints how many reads will be processed.
- Checks that tags are indeed made of ATGC bases.
- Added a check for cell barcode correction and collapsing threshold for given whitelist.
- New option to allow for a sliding window when aligning TAGS. `--sliding-window`
  Use when you have a variable sequence before your tags.

### Changed
- Sparse matrix is now based on int32 instead of int16. Fixing issue #40
- Cell barcode correction without a whitelist is not outputing any error anymore.
  This is due to the fact that it uses a hard threshold based on the `-cells` option
  if it kind find one. Fixing issues #29, #36
- Upgraded umi_tools to `1.0.0`
- fixed the error linked to umi_tools version update `getCellWhitelist` #45.


## [1.4.0] - 24.01.2019
### Added
- Enabled parallelization using multiprocessing. You can choose how many
  cores/threads you want to use with the `-T` `--threads` option. Takes max
  cpu by default. It will run on only one core if you run 1'000'000 reads or less.
- The output is now given in a gzipped mtx format. This is to ensure smooth usage for 
  new/larger datasets such as data from novaseq sequencers. For those who still want
  to use the dense format, there is an option `--dense` which will add the dense output
  as a tsv format.
  The sparse outputs are compatible the `Read10x` from the [Seurat](https://satijalab.org/seurat/) package.
- You get both umi and read counts as output. This can help with overamplification
  estimation.
- UMI and cell barcodes are now corrected based on [umi_tools](https://github.com/CGATOxford/UMI-tools).
- A small report is now produced as `report.yaml` giving some statistics about the run
  as well as which parameters have been used.
- Unittesting has been added using pytest. This should help prevent further
  unforseen bugs.
- New option for trimming the start of read2 sequences. `--start-trim`
- `--version` option now prints the version.

### Changed
- CITE-seq-Count uses less memory and is faster.
- Changed how you select the unmapped tags. The option is now `-ut` instead of `-uc`
  and gives you the top most unmapped tags.
- Expected cells, `-cells` is now required and not mutually exclusive with whitelists.
- Hamming distance option has been changed from `-hd` to `--max-error`.
  Please refer to the [documentation](https://hoohm.github.io/CITE-seq-Count/) for more information.

### Removed
- The default dense output matrix is gone and replaced by the sparse equivalent.
  See the [documentation](https://hoohm.github.io/CITE-seq-Count/) on how to read it.
- Regex is gone. Finding out what doesn't map is given by the unmapped file option.


## [1.3.4]

### Added
- Stripping of numbers and dash (-) to the *whitelist* barcodes so the 10X
  `barcodes.tsv` file could be directly used.
- When TAGs are too close based on the maximum Levenshtein distance allowed,
  print the offending pairs with its distance to identify them.
- When storing unmatched tags, added the possibility to specify a minimum
  counts cutoff to avoid printing every unmatched possibility (option -uc).
- Documentation to the functions.

### Changed
- `Levenshtein.hamming` usage was replaced by `Levenshtein.distance` in order
  to keep the sequences with INDELs matching against the TAGs.
- The redundant classifications `no_match` and `bad_struct` were merged into one
  (`no_match`); now, sequences match or do not.
- Decreased the running time by compiling the regex.
- File handles are released once done reading Read1/Read2 files.
- The different length TAGs are now being matched: first, by the regex pattern,
  which will always pick the longest match within the maximum distance allowed;
  second, by looping through the TAGs sorted by decreasing length, and choosing
  the first match within the maximum distance allowed.
- General code improvements and optimisations: releasing file handles sooner,
  compartimentalising code, etc.

### Removed
- Removed `ambiguous` classification; it's not a possibility because it's being
  handled when checking the distance between TAGs. E.g.
  Lets assume `max_hamming_distance = 2`, then:
  if any two TAGs are 2 hamming distance away from each other, then the program
  aborts with a message; thus, we could never be finding two sequences being
  two hamming distance away from more than one TAG.

### Fixed
- Reduced the memory usage by 1/3 by removing the `UMI_reduce` set; because it
  was based on the `unique_lines` set, where Read1 is already trimmed to
  Barcode+UMI length, it was providing no gain. (Issue #17).
- Hamming distance was not being properly applied because the equal sign was
  missing from the the conditional (`if min_value >= args.hamming_thresh`).
- When applying the regular expression to filter patterns to keep, by using the
  `i` in the fuzzy logic it was only accepting `insertions` as mismatches, not
  allowing proper mismatches or deletions. The `i` was changed to `s`, which
  means `substitution` in general.
- The `merge` technique for creating the common regex pattern for all the
  requested ABs/TAGs was replaced by a simple `in list` technique; the `merge`
  technique was faulty, matching sequences that are clearly different
  (`bad_struct`), being later classified as `no_match` instead. For example:
  Tag1      AAAAAA
  Tag2      TTTTTT
  Pattern   [AT][AT][AT][AT][AT][AT]
  Read2     ATATAT...
  In this scenario, Read2 is being matched and classified as `no_match` instead
  of `bad_struct`.


## [1.3.2] - 2018-10-22
### Added
- Printing version now from the help
- Added a possibility to store unmatched tags in a file using the option `-u`


## [1.3] - 2018-09-05
### Added
- `-l` `--legacy` option now available. This option will then end of the regex
  to find the TAG barcode in R2. Use `-l` if you have TAG barcodes ending with
  [TGC] + polyA tails.
- New warning that checks R1 length in the fastq file and the cell and umi
  barcodes given as input.

### Changed
- TAG barcode sequence are now added to the name of the tag in the rows of the
  results. This will help reproducibility since the barcode sequence will
  already be in the count results.


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