---
layout: default
title: PepSIRF Tutorials
permalink: /pepsirf-tutorials/
---
### PepSeq Manual Outline

Demultiplex the high-throughput sequencing data and assign raw read counts for each peptide. **[PepSIRF demux]**.

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
--sindex1 Index1 \
--sindex2 Index2 \
-d diagnostics.out
```

Normalize read counts (reads per million, RPM) in the demultiplexed data to account for variability in sequencing depth between samples. **[PepSIRF norm, col_sum]**.

```
pepsirf norm -a col_sum -p comboWR_raw_2mm_i1mm.tsv -o comboWR_raw_2mm_i1mm_CS.tsv >> norm.out
```

Select negative control samples and generate groups of peptides with similar abundances to be used for  Z score calculations. **[PepSIRF bin]**.

&ensp;&ensp;&ensp;&ensp;Create a .txt file with the names of the negative control samples. One sample name per line (*neg_control_names.txt*).

&ensp;&ensp;&ensp;&ensp;Use pepsirf subjoin to select negative control samples identified in neg_control_names.txt from column sum normalized dataset.

```
pepsirf subjoin -i comboWR_raw_2mm_i1mm_CS.tsv,neg_control_names.txt -o comboWR_raw_2mm_i1mm_CS_neg_ctrl_only.tsv >> subjoin.out
```

&ensp;&ensp;&ensp;&ensp;Use column sum normalized negative control matrix to create bins.

```
pepsirf bin -s comboWR_raw_2mm_i1mm_CS_neg_ctrl_only.tsv -b 300 -r 1 -o comboWR_raw_2mm_i1mm_b300r1_bins.tsv >> bin.out
```

Further normalize RPM counts from step 128 by subtracting the average RPM across negative controls. **[PepSIRF norm, diff]**.

```
pepsirf norm -a diff -p comboWR_raw_2mm_i1mm_CS.tsv -o comboWR_raw_2mm_i1mm_SBD.tsv --negative_id SB  >> norm.out
```

Calculate Z scores for each peptide using the bins created in step 129. **[PepSIRF zscore]**.

```
pepsirf zscore -s comboWR_raw_2mm_i1mm_SBD.tsv -o comboWR_raw_2mm_i1mm_Z-HDI75.tsv -n comboWR_raw_2mm_i1mm_Z-HDI75.nan -b comboWR_raw_2mm_i1mm_b300r1_bins.tsv -d 0.750000 >> zscore.out
```

Generate lists of enriched peptides for each sample based on thresholds for Z score and/or other metrics of interest. **[PepSIRF enrich]**.

You will need to first create a tsv file with paired sample names. First get all sample names using PepSIRF info module. Then create a tab delimited file with column for replicate 1 and a column for replicate 2 (see *comboWR_raw_2mm_i1mm_PN.tsv* as example).

```
info -i comboWR_raw_2mm_i1mm.tsv -s comboWR_raw_2mm_i1mm_SN.tsv >> info.out
```
<br>
```
pepsirf enrich -t comboWR_raw_2mm_i1mm_thresh.tsv -s comboWR_raw_2mm_i1mm_PN.tsv -r comboWR_raw_2mm_i1mm.tsv --raw_score_constraint 250000 -x _enriched.txt -f enrichFailReasons.tsv -o 10Z-HDI75_100CS_250000raw >> enrich.out
```
