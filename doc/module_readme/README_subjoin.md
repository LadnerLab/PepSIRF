#### Module Subjoin
```
PepSIRF: Peptide-based Serological Immune Response Framework subjoin module:
  -h [ --help ]                      Produce help message
                                     
  --names_list arg                   The name of the file containing the names 
                                     of the peptides to include, one per 
                                     line.Each of the peptides in this list 
                                     should be a row in the matrix provided in 
                                     the --scores argument. 
                                     
  --scores arg                       The name of the file containing scores for
                                     peptides. This should be a matrix peptide 
                                     name row labels and samplenames for column
                                     labels. Most probably this will be the 
                                     output generated by the normalize or demux
                                     modules.
                                     
  --output arg (=subjoin_output.tsv) The name of the file to write output 
                                     scores to. The output will be in the form 
                                     of the input, but with only the peptides 
                                     found in the namelists. 
```