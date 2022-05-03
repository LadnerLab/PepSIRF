---
layout: default
title: PepSIRF Tutorials
permalink: /pepsirf-tutorials/
---
### PepSeq Manual Outline

Demultiplex the high-throughput sequencing data and assign raw read counts for each peptide.  In this example read 1 and 2 are input as gunzipped fastq files. The file containing the barcodes used in the sequencing run are input to the index flag. The list of samples linked to the barcodes is input to the sample list flag. A fasta file containing the nt sequences of the peptides without the adaptor sequences is passed to the library flag. In this example the unique DNA tag sequence starts at the 43rd nt, is 90nt long, and we are allowing for 2 mismatches in the DNA tag sequence (\-\-seq 43,90,2). The same order of numbers is used to indicate the start, length and number of allowed mismatches for index1 and index2.

```
pepsirf demux --input_r1 sampled_R1.fastq.gz \
--input_r2 sampled_I1.fastq.gz \
--index BSC_FR_barcodes_Plus.fa \
-o demux_tutorial_raw_2mm_i1mm.tsv \
--samplelist ZZ_sample_list_PCV_Edit.tsv \
--library PCV_coded_hits.fna \
--read_per_loop 80000 \
--num_threads 12 \
--seq 43,90,2 \
--index1 12,12,1 \
--index2 0,8,1 \
-d diagnostics.out
```
<br>

Using the pepsirf norm module with the col_sum flag, normalize the demultiplexed read counts to reads per million (RPM) to account for variability in sequencing depth between samples. **[PepSIRF norm, col_sum]**.

```
pepsirf norm -a col_sum -p comboWR_raw_2mm_i1mm.tsv -o comboWR_raw_2mm_i1mm_CS.tsv >> norm.out
```
<br>

Select negative control samples from the full column sum normalized dataset and use them to generate groups of peptides (bins) with similar abundances to be used for  Z score calculations. **[PepSIRF bin]**.

+ Create a .txt file with the names of the negative control samples. One sample name per line (*neg_control_names.txt*)

+ Use pepsirf subjoin to select negative control samples identified in neg_control_names.txt from column sum normalized dataset

```
pepsirf subjoin -i comboWR_raw_2mm_i1mm_CS.tsv,neg_control_names.txt -o comboWR_raw_2mm_i1mm_CS_neg_ctrl_only.tsv >> subjoin.out
```
<br>

+ Use column sum normalized negative control matrix to create bins

```
pepsirf bin -s comboWR_raw_2mm_i1mm_CS_neg_ctrl_only.tsv -b 300 -r 1 -o comboWR_raw_2mm_i1mm_b300r1_bins.tsv >> bin.out
```
<br>

Further normalize the RPM counts by subtracting the average RPM across negative controls to account for variations in peptide abundance in the translated peptide library as well as any background binding of peptides to the capture beads. **[PepSIRF norm, diff]**.

```
pepsirf norm -a diff -p comboWR_raw_2mm_i1mm_CS.tsv -o comboWR_raw_2mm_i1mm_SBD.tsv --negative_id SB  >> norm.out
```
<br>

Calculate Z scores for each peptide using the negative control subtracted RPM matrix and the bins created from the RPM normalized negative control samples. **[PepSIRF zscore]**.

```
pepsirf zscore -s comboWR_raw_2mm_i1mm_SBD.tsv -o comboWR_raw_2mm_i1mm_Z-HDI75.tsv -n comboWR_raw_2mm_i1mm_Z-HDI75.nan -b comboWR_raw_2mm_i1mm_b300r1_bins.tsv -d 0.750000 >> zscore.out
```
<br>

Generate lists of enriched peptides for each sample based on thresholds for Z score and/or other metrics of interest. **[PepSIRF enrich]**.

+ You will need to first create a tsv file with paired sample names. First get all sample names using PepSIRF info module. Then create a tab delimited file with column for replicate 1 and a column for replicate 2 (see *comboWR_raw_2mm_i1mm_PN.tsv* as example).

```
pepsirf info -i comboWR_raw_2mm_i1mm.tsv -s comboWR_raw_2mm_i1mm_SN.tsv >> info.out
```
<br>
```
pepsirf enrich -t comboWR_raw_2mm_i1mm_thresh.tsv -s comboWR_raw_2mm_i1mm_PN.tsv -r comboWR_raw_2mm_i1mm.tsv --raw_score_constraint 15000 -x _enriched.txt -f enrichFailReasons.tsv -o 10Z-HDI75_100CS_15000raw >> enrich.out
```
<br>
