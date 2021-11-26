---
layout: default
title: Changelog
permalink: /changelog/
---

# Changelog

## 1.4.0 | 2021-07-09

Version 1.4.0 adds multiple features and one bug fix for s_enrich, p_enrich, and link. CMakelists has been updated and a new module ‘enrich’ has been introduced.


### New Features:

<strong>Module added: enrich (Issue #114).</strong> The p_enrich module was altered to allow for flexibility in the number of replicates for each sample and renamed ‘enrich’. This new module can now provide the functionality of both s_enrich and p_enrich, and therefore, these two modules will no longer be available. Additionally, this module is able to handle >2 replicates.

<strong>Enrich: new optional output file (Issue #103).</strong> An optional flag (-f, --enrichment_failure_reason) is now available. If used, a .tsv file will be generated to document each sample for which an enriched peptide file was NOT generated, as well as the reason  why.

<strong>CMakelists: Big Sur support (Issue #117).</strong> ‘-Xpreprocessor’ has been added to the command setting CMake C++ flags in order to support compilation on Mac OS Big Sur.


### Bug Fixes:

<strong>Link: Issue #116.</strong> A vague and system-dependent error occurred when --protein_file sequence names were not found in the --meta file. Modifications have been made to properly handle this situation and provide a clear and consistent error message.

## 1.3.7 | 2021-06-28
Version 1.3.7 adds one feature and one bug fix to norm.

### Bug Fixes:

<strong>Norm: Issue #125.</strong> When using a separate negative control matrix (--negative_control) for diff, ratio or diffratio normalizations, previous versions assumed that the order of the peptides (rows) was identical to the order in the primary data matrix (-p, --peptide_scores). This has been changed to properly account for any order in both rows and columns.

<strong>Norm: Issue #104.</strong> The norm module help message for option (--peptide_score, -p) has been updated.

## 1.3.6 | 2021-06-09
Version 1.3.6 adds several features and fixes several issues in demux, zscore, and subjoin.


### New Features:

<strong>Demux: new warning (Issue #96).</strong> The module now includes a warning for the user when index names from the (--samplelist, -s) file are not included in the index fasta file (--index, -i).

<strong>Zscore: new error handling (Issue #97).</strong> The (--bins, -b) file is now verified to be of the correct type. If the file is not of the correct format, the user will now receive an easily understandable error message.

### Bug Fixes:

<strong>Subjoin: Renaming bug (Issue #99).</strong> There was a bug related to the subjoin module renaming feature, which, under certain circumstances, resulted in the wrong samples being included in the output matrix. A fix has been made to prevent this feature from causing errors in subjoin.

<strong>Demux: Incorrect diagnostics calculations (Issue #100).</strong> The diagnostic features added in v1.3.2 miscalculated the number of index matches due to the module using every barcode contained in (-i, --index) for evaluation. Barcodes in the index file are now only used in the analysis if they are present in the (-s, --samplelist) file.


## 1.3.5 | 2021-03-02

<strong>Norm: Incorrect handling of input and output data (PR #93).</strong> This version changes how the negative control arguments are checked so that one of these flags (--negative_names) or (--negative_id) is only required when using the diff, ratio, or diff-ratio approaches.

<strong>Norm: Added missing transposing and updated approach using negative control matrix (PR #93).</strong> In v1.3.4, the optional negative control matrix (--negative_control) was not transposed before being accessed in calculations. By not including this, segfaults or incorrect averages were likely to occur.

<strong>P_enrich: Raw score threshold causing segfault (PR #94).</strong> In v1.3.4, the p_enrich module does not properly initialize the raw score threshold. The raw score thresholds no longer default to a container of 2 values both set to '0.0'. Raw score threshold will only contain the provided values through the (--raw_score_constraint) option for p_enrich and (--threshhold_file, -t) option.


## 1.3.4 | 2020-12-28

### New Features:

<strong>Demux: new sample list format (Issue #89).</strong> The demux module’s (--samplelist, -s) file must now include a header row and may optionally include more than the required 3 columns used for demultiplexing. In v1.3.4, the previous, header-less (--samplelist, -s) file format has been deprecated. The columns needed for demultiplexing are now identified by header names and therefore, order of these columns are no longer important. Header names can be specified using 3 new arguments:
 (--sname), (--sindex1) and (--sindex2).

<strong>Norm: new control-based approaches (Issue # 82).</strong> Three new normalization approaches were added( diff, ratio, and diff-ratio), all of which utilize column-sum normalized negative control values. To use these methods, the user must provide one of these two new arguments: (-s, --negative_id) or (-n, --negative_names). . If both are provided, then the (-n, --negative_names) will be used. A new optional argument (--negative_control) was added to allow the user to specify an independent data matrix containing the negative controls. 

<strong>P_enrich: increased threshold flexibility (Issue #84).</strong> Modified the way peptide-level thresholds are provided. Dedicated command line flags have been deprecated (--zscores, --zscore_constraint, --norm_scores, --norm_score_constraint). Data files and threshold values are now all specified using a single tab-delimited file (--threshold_file), which has the following format per row: a data matrix filename and either a single threshold or a comma-delimited pair of thresholds. This file should contain one row per data matrix, and any number of data matrices can be included.

## 1.3.3 | 2020-10-16

### New Features:

<strong>Subjoin: Adding second option for providing input data (PR #86).</strong> Two separate options have replaced the single filter scores option as input approaches. (--multi_file, -m) is for providing multiple input data matrices and name lists in a single tab-delimited file, while (--input, -i) is for providing a single data matrix and name list pair directly on the command line.

### Bug Fixes:

<strong>Subjoin: Remove "Sequence Names" from input matrix datasets(Issue #83).</strong> Subjoin no longer reads the first column of the first row when validating the sample names in the header row.


## 1.3.2 | 2020-09-20
Version 1.3.2 added features to demux and zscore as well as fixed one bug with compilation of a package for PepSIRF.

### New Features:

<strong>Demux: Additional diagnostic information (Issue #34).</strong> The demux module now features additional output to aid in output diagnostics. The percentage of index 1 matches, index 2 matches, and the variable region matches are output at the end of the run to the terminal.

<strong>Demux: Additional option for diagnostic output (Issue #34).</strong> The (--diagnostic_info) option accepts a filename to output diagnostic information to. The output contains sample names in the first column, and the following columns contain the index1, index2 and variable region matches.

<strong>Zscore: Additional feature for calculating z score (Issue #78).</strong> A highest density interval approach is now a feature of the z_score module.

### Bug Fixes:

<strong>PepSIRF: Fixing compatibility issue between ZLIB and Mac OS (Issue #72).</strong> ZLIB is now functional on both Mac and Linux systems as of this update.


## 1.3.1 | 2020-07-09
Version 1.3.1 fixes an issue with s_enrich and disabled a nonfunctioning utility for Mac operating systems used in PepSIRF.

### Bug Fixes:

<strong>S_enrich: Fixing data element mismatches during analysis (Issue #73).</strong> The sample names in the output of s_enrich were being switched around. The sample names are no longer mismatched with their data.

<strong>PepSIRF: Disabling ZLIB feature for Mac users (Issue #72).</strong> When attempting to compile PepSIRF with this ZLIB on Mac systems, a compilation error would occur. This feature will be temporarily disabled for Mac users until a fix is found.


## 1.3.0 | 2020-06-22
Version 1.3.0 adds features to deconv, demux, info, link, p_enrich, and subjoin. The help info was updated for multiple modules. It also includes bug fixes to bin, s_enrich, and the PepSIRF testing executable.

### New Features:

<strong>Deconv: Standardized format of "--linked" file for input (Issue #64).</strong> The Deconv module now includes a required option (-l,--linked). This is a file name with the file formatted in the output provided by the Link module.

<strong>Deconv: Simplified scoring method selection (Issue #61).</strong> The Deconv module now uses a single scoring strategy option. The strategies will no longer be separate options; now the new option (--scoring_strategy) accepts one of three string inputs: ‘summation’, ‘integer’, or ‘fraction’. By default, ‘summation’ is used.

<strong>Demux: Providing error messages (Issue #30).</strong> When parsing samples or fasta files in demux, there will now be a runtime error thrown when the file fails to be opened properly.

<strong>Info: Increasing significant digits (Issue #51).</strong> The output from the info module is now in fixed point notation with set precision at 2 decimal places.

<strong>Link: Switch source of taxonomic info to metadata file (Issue #50).</strong> The link module now uses a single metadata file to obtain taxonomic information for use in the module. One column for the protein sequence names and one column for the metadata to be used in generating the linkage map. Multiple columns of metadata can exist in the map, but only one ID column can be used.

<strong>P_enrich: Differentiate two "-s" flags (Issue #55).</strong> The single hyphen flag (-s) for outfile suffix option (--outfile_suffix) and the option (--samples) have been disambiguated by distinguishing the (-s) flag (--samples, -s) and (--outfile_suffix, -x).

<strong>Subjoin: Allow subjoin without a namelist (Issue #26).</strong> The subjoin module no longer requires a name list of sample names to be included with the score matrix for the (--filter_scores) option.

<strong>PepSIRF: Update to help info (Issue #53).</strong> The help info for the bin, deconv, demux, info, link, norm, p_enrich, s_enrich, subjoin, and z_score modules provided by (-h) have all been refactored to improve readability and accuracy.

<strong>PepSIRF: "pepsirf_test" executable does not compile (Issue #47).</strong> Unit tests in the testing build have been updated to reflect the changes made in the various modules for this update.

### Bug Fixes:

<strong>Bin: Last bin is smaller than specified "--bin_size" (Issue #58).</strong> When there is not a large number of peptides with zero count, the last bin can result in fewer than the number of peptides specified with the (--bin_size) option. The bin size increases by one until the smallest bin is equal to the size of the minimum bin size.

<strong>S_enrich: Minimum z score not properly set (Issue #54).</strong> The s_enrich module in previous versions was sharing the same minimum threshold provided for zscore for both the norm and zscore thresholds. This has been fixed so the thresholds are correctly used based on what is provided by the min score options.

