#### Module Subjoin
```
PepSIRF: Peptide-based Serological Immune Response Framework subjoin module:
  -h [ --help ]                      Produce help message
                                     
  --filter_scores arg                Comma-separated filenames (For example: 
                                     score_matrix.tsv,peptide_names.txt ) for a
                                     score matrix and a file containing the 
                                     names of peptides to keep in the score 
                                     matrix. The score matrix should be of the 
                                     format output by the demux module, with 
                                     sample names on the columns and peptide 
                                     names on the rows. The peptide namelist 
                                     must have one name per line. To use 
                                     multiple name lists with multiple score 
                                     matrices, include this argument multiple 
                                     times.
                                     
  --output arg (=subjoin_output.tsv) The name of the file to write output 
                                     scores to. The output will be in the form 
                                     of the input, but with only the peptides 
                                     found in the namelists. 
                                     
```
