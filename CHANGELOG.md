# Change Log
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/)
and this project adheres to [Semantic Versioning](http://semver.org/).

## [1.4.0]
### Added
- Enabled parallelization using multiprocessing. You can choose how many
  cores/threads you want to use with the `-T` `--threads` option. It takes max
  cpu by default.
- The output is now given in a gzipped mtx format. This is to ensure smooth usage for 
  new/larger datasets such as data from novaseq sequencers.
  It is also compatible with the `Seurat` package `Read10x` function.
- You get both umi and read counts as output. This can help with overamplification
  estimation.
- UMI and cell barcodes are now corrected based on umi_tools.
- A small report is now produced as `report.yaml` giving some statistics about the run
  as well as which parameters have been used.
- Unittesting have been added now using pytest. This should help prevent further
  unforseen bugs.

### Changed
- CITE-seq-Count uses less memory now.
- Changed how you select the unmapped tags. The option is now `-ut` instead of `-uc`
  and gives you the top most unmapped tags.
- Expected cells, `-cells` is now required and not mutually exclusive with whitelists.
- Hamming distance option has been changed from `-hd` to `--max-error`.
  Please refer to the documentation for more information.

### Removed
- Dense output matrix is gone and replaced by the sparse equivalent.
- pandas dependency.
- Regex are gone. Finding out what doesn't map is given by the unmapped file option.


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
  allowing proper mismatches or deletions. The `i` was changed to `e`, which
  means `error` in general.
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