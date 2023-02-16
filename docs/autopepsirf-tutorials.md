---
layout: default
title: AutoPepSIRF Tutorials
permalink: /autopepsirf-tutorials/
---
### PepSeq AutoPepSIRF Outline

Demultiplex the high-throughput sequencing data and assign raw read counts for each peptide.  In this example read 1 and 2 are input as gunzipped fastq files. The file containing the barcodes used in the sequencing run are input to the index flag. The list of samples linked to the barcodes is input to the sample list flag. A fasta file containing the nt sequences of the peptides without the adaptor sequences is passed to the library flag. In this example the unique DNA tag sequence starts at the 43rd nt, is 90nt long, and we are allowing for 2 mismatches in the DNA tag sequence (\-\-seq 43,90,2). The same order of numbers is used to indicate the start, length and number of allowed mismatches for index1 and index2. The (\-\-sindex) flags tell the program the names of the column headers in the samplelist file that indicate the identity of the indexes associated with each sample. **[PepSIRF demux]**.

In this example, our [input directory](https://github.com/LadnerLab/PepSIRF/tree/master/tutorial_files/demux) will contain all of the files required to demultiplex the high-throughput sequencing data and assign raw read counts.

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

For the remainder of the tutorial the [input directory](https://github.com/LadnerLab/PepSIRF/tree/master/tutorial_files/autopepsirf) will contain the raw read files needed. Examples of the expected outputs can be found in the [expected output directory](https://github.com/LadnerLab/PepSIRF/tree/master/tutorial_files/autopepsirf/expected_outputs)

Using the pepsirf norm module with the col_sum flag, normalize the demultiplexed read counts to reads per million (RPM) to account for variability in sequencing depth between samples. **[PepSIRF norm, col_sum]**.

```
pepsirf norm -a col_sum -p comboWR_raw_2mm_i1mm.tsv -o comboWR_raw_2mm_i1mm_CS.tsv >> norm.out
```
<br>

Select negative control samples from the full column sum normalized dataset and use them to generate groups of peptides (bins) with similar abundances to be used for  Z score calculations. **[PepSIRF bin]**.

+ Create a .txt file with the names of the negative control samples. One sample name per line (*neg_control_names.txt*)

+ Use pepsirf subjoin to select negative control samples identified in *neg_control_names.txt* from column sum normalized dataset

```
pepsirf subjoin -i comboWR_raw_2mm_i1mm_CS.tsv,neg_control_names.txt -o comboWR_raw_2mm_i1mm_CS_neg_ctrl_only.tsv >> subjoin.out
```
<br>

+ Use column sum normalized negative control matrix to create bins

```
pepsirf bin -s comboWR_raw_2mm_i1mm_CS_neg_ctrl_only.tsv -b 300 -r 1 -o comboWR_raw_2mm_i1mm_b300r1_bins.tsv >> bin.out
```
<br>

Run autopepsirf Qiime2 plugin to automatically run PepSIRF normalization, pair names file creation, z score calculation, and peptide enrichment analysis from one command. Click [here](https://ladnerlab.github.io/pepsirf-q2-plugin-docs/Pepsirf_Plugins/q2-autopepsirf/#installation) for instructions on the installation of the autopepsirf Qiime2 plugin. Here we are using the (\-\-p\-negative\-id) flag to show that any sample name starting with “SB” is to be considered a negative control sample. For the calculation of enriched peptides the Z score threshold has been set to 10 and HDI to 75%. The (\-\-p\-infer\-pairs\-source) flag will look for sample pairs that match except for after the last underscore, so sample1_A and sample1_B would be treated as a replicate pair.

```
qiime autopepsirf diffEnrich-tsv \
--i-raw-data comboWR_raw_2mm_i1mm.tsv \
--p-negative-id SB \
--i-bins comboWR_raw_2mm_i1mm_b300r1_bins.tsv \
--p-raw-constraint 15000 \
--p-exact-z-thresh 10 \
--p-hdi 0.75 \
--p-infer-pairs-source \
--p-pepsirf-binary pepsirf \
--p-pepsirf-tsv-dir ./aps-TSV \
--p-tsv-base-str comboWR_raw_2mm_i1mm \
--output-dir autopepsirf-diffEnrich 
```
<br>