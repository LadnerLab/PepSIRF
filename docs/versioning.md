---
layout: default
title: See What's New
permalink: /changelog/
---

# What's New: 

### [Version 1.3.5 (03.02.2021)](#135)
### [Version 1.3.4 (12.28.2020)](#134)
### [Version 1.3.3 (10.16.2020)](#133)
### [Version 1.3.2 (09.20.2020)](#132)
### [Version 1.3.1 (07.09.2020)](#131)
### [Version 1.3.0 (06.22.2020)](#130)

## 1.3.5
<strong>Overview:</strong> Version 1.3.5 fixes several issues found in the norm and p_enrich modules.
<details><summary>Details</summary>
    <br>
    <strong>Norm bug fixes.</strong> In version 1.3.4, three new approaches were added: diff, ratio, and diff-ratio. To use these normalization approaches, the user must provide either a negative control name list (--negative_names) or a negative control ID string (--negative_id). Therefore, a check was added in v1.3.4 to ensure that one of these was provided; however, this check was implemented for all normalization approaches, even those that do not utilize negative controls. v1.3.5 changes this behavior so that one of these flags is only required when using the diff, ratio, or diff-ratio approaches. Additionally, in v1.3.4 the optional negative control matrix (--negative_control) was not transposed before being accessed in calculations. By not including this, segfaults or incorrect averages were likely to occur. This has been fixed in v1.3.5.<br>
    <br>
    <strong>P_enrich bug fixes.</strong> In v1.3.4, the p_enrich module did not properly initialize the raw score threshold or each matrix file threshold. The threshold contained by default two empty zero values where, if used, the raw score constraint would overwrite one zero, but leave the other. The same would happen for every input matrix when only a single threshold was provided. This resulted in incorrect threshold comparisons with the sample pairs. The raw score constraint and matrix thresholds now only contain the provided values through their respective p_enrich raw score constraint option and threshold file filename,threshold pairs.<br>
    <br>

</details>

## 1.3.4
<strong>Overview:</strong> Version 1.3.4 adds new features to p_enrich, norm, and demux.
<details><summary>Details</summary>
    <br>
    <strong>Demux feature added.</strong> The demux module’s sample list file may now include more than 3 columns for demultiplexing. Prior to v1.3.4 the samplelist has been restricted to a particular tab-delimited, ordered format with the first column as an index1 column, the second as an optional index2 column, and the last column as a sample name column. In v1.3.4, this format has been deprecated for a more flexible and scalable replacement. The sample list is still a tab-delimited file but now with a header row that utilizes 3 new flags. The header row of the sample list must contain at least one sample name header and one index name header for identification of the columns that will be used in the run. The header name that represents the samplename column is specified by a samplename flag (--sname). By default, this is set to 'SampleName'. There are two index header name flags that each represent an index column. The index one header name (--sindex1) and the optional index two header name (--sindex2). By default, index one is set as 'Index1' and index two is set as 'Index2'. There is no restriction on the number of columns in the sample list, but every column must have a header name. References to forward and reverse index were also removed in comments throughout the code, but to no effect on the program from the user end.<br>
    <br>
    <strong>Norm feature added.</strong> The norm module has been updated to include three new approaches: diff, ratio, and diff-ratio. For these three approaches, a negative control average is calculated by summing the scores of each valid sample of each peptide and dividing by the total samples of that associated peptide. To determine valid samples for the peptides in the matrix, a filtering approach is required as a norm option. Currently available is the negative id or negative names list where only one or the other can be used, but if both are provided, then the negative names list is used over the id. The negative id is a unique substring prefixed to the sample names in the matrix. The program will look through the sample names and select any sample name that starts with that substring. The negative names list is a comma-separated list of full sample names provided as the input for the norm option. All sample names from the list that are found in the matrix will be included. The diff approach normalizes the scores using the difference method. For each sample of each peptide, the difference of the (sample,peptide) score and the peptides negative score average is calculated. The ratio approach normalizes the scores using the ratio method. For each sample of each peptide, the ratio of the (sample,peptide) score and the peptides negative score average is calculated. The diff-ratio approach normalizes the scores using the combination of the difference and ratio methods. For each sample of each peptide, the difference of the (sample,peptide) score and the peptides negative score average is calculated and then the ratio of the difference and the peptides negative score average is calculated.<br>
    <br>
    <strong>p_enrich feature added.</strong> The p_enrich module now uses a tab-delimited file which has the following format: a matrix filename tab-delimited by a single or pair of thresholds, comma-delimited if a pair. The matrix file no longer requires separate input options to provide matrices as well as no longer holds a restriction on the number of matrices that can be provided.<br>
    <br>

</details>

## 1.3.3
<strong>Overview:</strong> Version 1.3.3 added some features and fixed an issue with the Subjoin module.
<details><summary>Details</summary>
    <br>
    <strong>Subjoin feature added.</strong> The filter scores option that takes two comma-delimited filenames or a single tab-delimited file containing comma-delimited files has been updated. To improve flexibility two separate options have replaced the single filter scores option as input approaches. The muli file option is one of these; this option accepts the name of a tab-delimited file containing a score matrix filename and a filename containing a list of sample names. The sample names list may have two tab-delimited columns. The second column names will replace the names of the sample names in the subjoin output. The second new option is input. This option now only accepts a single comma-delimited pair of filenames, the matrix score filename and the list of sample names filename. The multi file flag takes in the same types of input, but as a line in a file and with the allowance of multiple lines. Both of these flags can be used together and both of these flags do not require the sample list filename. If the filename is excluded, then all samples in the score matrix will be used.<br>
    <br>
    <strong>Subjoin bug fixed.</strong> Previously, when you did not provide a sample name list, the program would warn you that “Sequence name” was not found in the input matrix. This was due to an unnecessary verification of the first column in the sample name header row. Subjoin no longer reads the first column of the first row when validating the sample names in the header row.<br>
    <br>
    
</details>

## 1.3.2
<strong>Overview:</strong> Version 1.3.2 added features to Demux and Zscore as well as fixed one bug with compilation of a package for PepSIRF.
<details><summary>Details</summary>
    <br>
    <strong>Demux feature added.</strong> The Demux module now features additional output to aid in output diagnostics. The percentage of index 1 matches, index 2 matches, and the variable region matches are output at the end of the run to the terminal. A new, optional flag has been included to provide diagnostic information as well. The “--diagnostic_info” option accepts a filename to output diagnostic information to. The output is structured as 3-4 tab-delimited columns. The first column is sample names from input matrix, the second column is the number of index 1 matches, the third column is the number of index pair matches, dependent on the inclusion of index 2, and the fourth column is the number of index pair matches that also had a variable region match, dependent on the inclusion of a library file.<br>
    <br>
    <strong>Zscore feature added.</strong> A highest density interval approach is now a feature of the Zscore module. This is an alternative approach to the preexisting trim approach. This approach discards outliers prior to calculating the mean and standard deviation. If the HDI is provided with the trim approach, trim will be overridden. The HDI option accepts a decimal value to be used as a percentage highest density interval.<br>
    <br>
    <strong>PepSIRF bug fixed.</strong> An issue where the ZLIB package could not be compiled on Mac systems properly has been fixed in this version. This issue prevented gzipped file usage on Mac operating systems. ZLIB is used in PepSIRF to allow for compressed files to be provided as input. ZLIB is now functional on both Mac and Linux systems as of this update.<br>
    <br>

</details>

## 1.3.1
<strong>Overview:</strong> Version 1.3.1 fixes an issue with Senrich and disabled a nonfunctioning utility for Mac operating systems used in PepSIRF.
<details><summary>Details</summary>
    <br>
    <strong>Senrich bug fixed.</strong> The sample names in the output of Senrich were being switched around. The sample names are no longer mismatched with their data.<br>
    <strong>PepSirf bug fixed.</strong> The ZLIB package is a utility used to handle compressed .zip files as input. When attempting to compile PepSIRF with this package on Mac systems, a compilation error would occur. This feature will be temporarily disabled for Mac users until a fix is found.<br>
    <br>

</details>

## 1.3.0
<strong>Overview:</strong> Version 1.3.0 adds features to Deconv, Demux, Info, Link, p_enrich, and Subjoin. The help info was updated for multiple modules. It also includes bug fixes to Bin, Senrich, and the PepSIRF testing executable.
<details><summary>Details</summary>
    <br>
    <strong>Deconv feature added.</strong> The Deconv module now includes a required option ‘-l,--linked’. This input is in the form of a file name with the file formatted in the output provided by the Link module meaning a tab-delimited pair per line where the first column is a peptide name and the second is a Id and score delimited by a ‘:’.<br>
    <br>
    <strong>Deconv feature added.</strong> The Deconv module now uses a single scoring strategy option. The strategies will no longer be separate options; now the new option ‘--scoring_strategy’ accepts one of three string inputs: ‘summation’, ‘integer’, or ‘fraction’. By default, ‘summation’ is used.<br>
    <br>
    <strong>Demux feature added.</strong> When parsing samples in Demux, there will now be a runtime error thrown when the file fails to be open properly. Parsing fasta files will now also throw a runtime error when the file cannot be opened properly.<br>
    <br>
    <strong>feature added.</strong> The output from the Info module is now in fixed point notation with set precision at 2 decimal places.<br>
    <br>
    <strong>Link feature added.</strong> To increase flexibility and simplicity, the Link module now uses a single metadata file to obtain taxonomic information for use in the module. The file includes one column for the protein sequence names and one column for the metadata to be used in generating the linkage map. Multiple columns of metadata can exist in the map, but only one ID column can be used. The input provided in the option is formatted as 3 comma-delimited values: the filename for the metadata file, the header for the column containing the protein sequence names, and the header for the column containing the ID to be used in creating the linkage map.<br>
    <br>
    <strong>p_enrich feature added.</strong> The single hyphen flag ‘-s’ for outfile suffix option ‘--outfile_suffix’ and the option ‘--samples’ have been disambiguated by distinguishing the ‘-s’ flag. The ‘-s’ flag now refers to the ‘--samples’ option while the new flag ‘-x’ will now be used with the ‘--outfile_suffix’ option.<br>
    <br>
    <strong>Subjoin feature added.</strong> The Subjoin module no longer requires a name list of sample names to be included with the score matrix for the ‘--filter_scores’ option. By excluding the name list filename in a tab-delimited input file or comma-delimited pair, every sample name will be included in the output.<br>
    <br>
    <strong>PepSIRF documentation update.</strong> The help info for the Bin, Deconv, Demux, Info, Link, norm, p_enrich, Senrich, Subjoin, and Zscore modules provided by ‘-h’ have all been refactored to improve readability and accuracy.<br>
    <br>
    <strong>Bin bug fixed.</strong> In cases in which there is not a large number of peptides with zero count, the last bin can result in fewer than the number of peptides specified with the ‘--bin_size’ option. This has been fixed by increasing the minimum bin size by one until there is no last bin with leftover peptides.<br>
    <br>
    <strong>Senrich bug fixed.</strong> The Senrich module in previous versions was sharing the same minimum threshold provided for zscore for both the norm and zscore thresholds. This has been fixed so the thresholds are correctly used based on what is provided by the min score options.<br>
    <br>
    <strong>PepSIRF testing update.</strong> Unit tests in the testing build have been updated to reflect the changes made in the various modules for this update. Input options, containers, variables, and function arguments required updates to continue functional testing of all PepSIRF modules.<br>
    <br>
</details>
