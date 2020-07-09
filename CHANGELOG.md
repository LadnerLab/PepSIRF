# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]
- Fixed bug causing sample names to be mismatched for s_enrich module output files.

## [1.3.0] - 2020-6-22
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
