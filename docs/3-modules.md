---
layout: default
title: Modules
permalink: /modules/
---

## Modules

### Bin
The bin module is used to create groups of peptides with similar starting abundances (i.e. bins), based on the normalized read counts for >= 1 negative controls. These bins can be provided as an input for zscore calculations using the zscore module.

### Deconvolution
The deconv module converts a list of enriched peptides into a parsimony-based list of likely taxa to which the assayed individual has likely been exposed. This module has two modes: batch and singular. In batch mode, the input given to '--enriched' is a directory containing files of enriched peptides for each sample (e.g., as output by _enrich or p_enrich). In this case, '--output' is a directory where the output deconvolution reports will be written. In singular mode, both '--enriched' and '--output' are treated as files, not directories. The chosen mode is determined by the type of argument provided with '--enriched'. If a directory is specified, batch mode will be used. If a file is specified, singular mode will be used.

### Demultiplexing
This module takes the following parameters and outputs counts for each reference sequence (i.e. probe/peptide) for each sample. For this module we define 'distance' as the Hamming distance D between a reference sequence r and a read sequence s. If D( r, s ) <= max_mismatches we say that s maps to r. Note that if multiple reference sequences, r1 and r2, are similar to s, D( r1, s ) <= max_mismatches and D( r2, s ) <= max_mismatches, then we say that s maps to the reference with the minimum distance. Additionally if D( r1, s ) == D( r2, q ) then we discard s as we cannot say whether s maps to r1 or r2.

### Info
This module is used to gather information about a score matrix. By default, the number of samples and peptides in the matrix will be output. Additional flags may be provided to extract different types of information. Each of these flags should be accompanied by an output file name, to which the information be written.

### Link
The link module is used to create the "--linked" input file for the deconv module. The output file from this module defines linkages between taxonomic groups (or other groups of interest) and peptides based on shared kmers.

### Normalize
The norm module is used to normalize raw count data to allow for meaningful comparison among samples.

### Paired (Duplicate) Enrichment
The p_enrich module determines which peptides are enriched in samples that have been assayed in duplicate, as determined by user-specified thresholds. Thresholds are rovided as comma-delimited pairs. In order for a peptide to be considered enriched, both replicates must meet or exceed the lower threshold and at least one replicate ust meet or exceed the higher threshold, independent of order. Note that a peptide must meet each specified threshold (e.g., zscore, norm count and raw count) in order to be considered enriched.

### Single-Replicate Enrichment
The s_enrich module determines which peptides are enriched in each sample, as determined by user-specified thresholds. Note that a peptide must meet each specified hreshold (e.g., zscore, norm count and raw count) in order to be considered enriched. This module will batch process all samples contained within a matrix and will enerate one output file per sample.

### Subjoin
The subjoin module is used to manipulate matrix files. This module can create a subset of an existing matrix, can combine multiple matrices together or perform a combination of these two functions.

### Zscore
The zscore module is used to calculate Z scores for each peptide in each sample. These Z scores represent the number of standard deviations away from the mean, with the mean and standard deviation both calculated separately for each bin of peptides.
