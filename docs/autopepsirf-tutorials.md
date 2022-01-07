---
layout: default
title: AutoPepSIRF Tutorials
permalink: /autopepsirf-tutorials/
---
### PepSeq AutoPepSIRF Outline

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

&ensp;&ensp;&ensp;&ensp;Use pepsirf subjoin to select negative control samples identified in *neg_control_names.txt* from column sum normalized dataset.

```
pepsirf subjoin -i comboWR_raw_2mm_i1mm_CS.tsv,neg_control_names.txt -o comboWR_raw_2mm_i1mm_CS_neg_ctrl_only.tsv >> subjoin.out
```

&ensp;&ensp;&ensp;&ensp;Use column sum normalized negative control matrix to create bins.

```
pepsirf bin -s comboWR_raw_2mm_i1mm_CS_neg_ctrl_only.tsv -b 300 -r 1 -o comboWR_raw_2mm_i1mm_b300r1_bins.tsv >> bin.out
```

Run autopepsirf.py wrapper to automatically run PepSIRF normalization, pair names file creation, z score calculation, and peptide enrichment analysis from one command.

```
/PepSIRF/extensions/autoPepSIRF.py \
-r comboWR_raw_2mm_i1mm.tsv \
--negative_id SB \
-b comboWR_raw_2mm_i1mm_b300r1_bins.tsv \
--rawThresh 250000 \
--zThresh 10 \
--hdi 0.75 \
--inferPairs \
>>pepsirf.out
```