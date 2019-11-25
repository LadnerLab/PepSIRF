#### Module Demux
```

PepSIRF: Peptide-based Serological Immune Response Framework score normalization module. 
:
  -h [ --help ]                         Produce help message
                                        
  -t [ --threads ] arg (=2)             The number of threads to use when 
                                        performing analyses.
                                        
  --peptide_scores arg                  Name of file containing peptide scores.
                                        This file should be tab-delimited with 
                                        the first column being peptide names, 
                                        and every next column should be 
                                        the peptide's score within a given 
                                        sample (the first item in the column). 
                                        This is exactly the format output by 
                                        the deconv module.
                                        
  -o [ --output ] arg (=norm_output.tsv)
                                        The name of the file to write output 
                                        to. The output is formatted in the same
                                        way the input 'peptide_scores' are 
                                        formatted, i.e. a score matrix with 
                                        samples on the columns and scores for a
                                        certain peptide on the rows. The score 
                                        for each peptide in the output has been
                                        normalized in the manner specified.
```
