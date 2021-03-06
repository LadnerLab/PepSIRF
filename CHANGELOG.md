# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [1.4.0] - 2021-07-09
- #117, CMakelists has been updated to include a new flag for the CXX flags: '-Xpreprocessor'. This flag is used to make compilation in different environments for cpp easier. This issue arose when pepsirf was attempted to be compiled in 'Big Sur' and failed to compile due to an error with '-fopenmp'.
- #116, link module had occuring error when protein sequences were not found in the metadata map. This has been changed so the situation is handled and an error is thrown stating a sequence was not found in the metadata file.
- #114, s_enrich and p_enrich have been merged into a single module 'enrich'. A single and pair of samplenames work with the same behavior as s_ and p_enrich respectively. Additionally, >2 replicates can be analyzed to generate enriched peptides. See help options for an update on the options.
- #103, enrich (previously two separate modules, s_enrich and p_enrich) module now features an optional flag that outputs in a tsv failed enrichment replicate sets. Any set replicates which failed to generate enriched peptides is listed in this file. Each row contains a column for replicate samplenames and a column containing the reason.

## [1.3.7] - 2021-06-28
- #125, norm module incorrectly stored peptide names with the assumption that in diff, diff-ratio, or ratio the control and original matrix are in the same order. The peptide names should not be assumed to be in the same order. This has been updated by changing the type of container used and the method of access.
- #104, norm module help message updated. ('-p', '--peptide_scores') option in the norm help message states in the final sentence "This file should be in the same format as the output from the deconv module.". This is incorrect - deconv should be demux - it should read "This file should be in the same format as the output from the demux module.".

## [1.3.6] - 2021-06-09
- #96, demux module now includes a warning when index/barcode names from the samplelist are not included in the fasta file provided by (--index). The warning includes a list of missing names.
- #97, zscore module now verifies the correct type of file is provided for (--bins). Assuming the incorrect file provided will be a score matrix tsv, the verification process is a check of the second line in the tsv. If a numerical value or 'inf' or 'NaN' is found, then an error is thrown with a message stating to check the file provided.
- #99, subjoin module bug fixed with names in output score matrix. The subjoin module features the ability to update the sample/peptide names for the output matrix using a second column in the namelist file provided as the second file in the pair. In version 1.3.5, this feature did not work correctly where name updates to existing names in the matrix that themselves would be updated, was prone to error by mixing up the associated score column to the name. eg. an original name "pv1_001" changed to "pv1_101" and and "pv1_001" being used elsewhere in the renaming was prone to switching names around due to stochastic ordering of names in the matrix. This has been fixed by adding an 'upate_labels' method and altered structure of the container storing names from the namelist file. There is also now a check that the namelist file can be opened before continuing the run.
- #100, demux module diagnostic feature added previously miscalculated counts for index1, pair, and var region matches. This was in part due to conditional checks for pair matches and the presence of barcode/index names from the index file that are not included in the samplelist file. The fix is to remove the unused names from the index file - only the generated container holding them, the file itself is not edited - and a more concise condition to verify a pair match occurred.

## [1.3.5] - 2021-03-02
- Fixed penrich bug with threshold verification for enrichment candidates and raw score pairs.
- Fixed bug with norm which forced input of id or names for use with negative control.
- Fixed bug with norm incorrectly accessing matrix elements for diff, ratio, diffratio.

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
