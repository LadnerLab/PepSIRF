# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]
- Norm no longer forcing input of id or names for use with negative control.

## [1.3.4] - 2020-12-28
- Penrich now uses a tab-delimited file containing a matrix file name and its threshold(s) per each line. A matrix file may contain zscores or normalized counts of each peptide.
- Demux sample list may now include more than 3 columns for demultiplexing. A header name must now be specified for each column in the list to specify the sample name column, index 1 column and potentially index 2 column. Updated comments to drop reference to forward and reverse indexes - now index 1 and index 2.
- Norm includes 3 new approaches: diff, ratio, and diff-ratio.

## [1.3.3] - 2020-10-16
- Subjoin bugfix where first line of list of provided matrix filenames may be skipped.
- Subjoin -f flag has become two separate flags, -i(--input) and -m(--multi_file). This is to increase flexability in
  providing a variety of input.
- Subjoin gave warnings when reading the "Sequence name" header as a sample name. This issue was created from 1.3.0 feature addition.

## [1.3.2] - 2020-09-20
- Fixed issue with Zlib compilation error on Mac.
- Added additional output features to demux to aid in diagnosis and tracing of the 2 possible indexes given and DNA tag matches.
- Added highest density interval as filtering option for zscore module.

## [1.3.1] - 2020-07-09
- Fixed bug causing sample names to be mismatched for s_enrich module output files.
- ZLib disabled for Mac OS temporarily to avoid compilation bug.

## [1.3.0] - 2020-06-22
- Deconv module now requires --linked file to be in format provided by link module output file.
- Deconv Module now uses a single scoring strategy flag that takes the name of the strategy as an argument.
- Updated help info provided by modules -h flag.
- Fixed bug where bin module last bin size falls below minimum.
- Changed link modules scoring strategy selection process to now use only one flag to which the strategy must be entered.
- Added functionality for metadata file usage in link module over an ID index number.
- Fixed bug where p_enrich -s arg referred to both the output suffix and the sample.
- Altered the info module output for sum of column scores to fixed-point notation.
- Fixed bug causing the max of the specified norm/score to be used in s_enrich
- Added option to subjoin where exclusion of a name list outputs all columns of a given matrix file.
- Added error checking and reporting to demux sample list parser and fasta parser.
- Fixed bug that prevented the 'pepsirf_test' executable from building.

## [1.2.2] - 2020-05-01
- (General) Added instructions for contributing to the PepSIRF software.

## [1.2.1] - 2020-04-06
- (Demux) Fixed a bug that allowed reads who had a forward index match but no reverse index match be output with a score of zero.
- (Demux) Fixed a bug causing a difference in scores between aggregate and translation-based non-aggregate scores.
## [1.2.0] - 2020-03-30
- Updated the help text of each module to display the current version number.
- Fixed demux flag names that were inconsistent between the CLI input and standard output.
- Demux can now read from gzipped fastq files
- Added a precision argument to norm that enables the specification of numeric precision in the output.
- Fixed bug caused by an incomplete Codon -> AA translation map

## [1.1.0] - 2020-02-26
- 2020-02-25: [General] Module help information now fills the entire line width on MacOS.
- 2020-02-25: [Demux] Added support for translation-based count aggregation.
- 2020-02-24: [Demux] Added support for reference-independent demultiplexing.
## [1.0.0] - 2020-02-14
Release of the first version of PepSIRF.
