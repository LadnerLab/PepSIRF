# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]
- Fixed bug that prevented the 'pepsirf_test' executable from building.

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
