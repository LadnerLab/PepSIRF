#### Module Demux
```
PepSIRF: Peptide-based Serological Immune Response Framework score normalization module. 
:
  -h [ --help ]                          Produce help message
                                         
  -p [ --peptide_scores ] arg            Name of file containing peptide scores. This file should be tab-delimited with the first column being peptide names, and every next column should be
                                          
                                         the peptide's score within a given sample (the first item in the column). This is exactly the format output by the deconv module.
                                         
  -c [ --col_sum ]                       Normalize the counts using a column-sum method. Output is the number of reads a peptide per reads million mapped. Note that if size_factors is also 
                                         included, the value of this flag will be ignored and the size_factors method is used. By default, col_sum normalization is used.
                                         
  -s [ --size_factors ]                  Normalize the counts using the size factors method (Anders and Huber 2010). Note that if this flag is included, the value of col_sum will be 
                                         ignored.
                                         
  --precision arg (=2)                   Output score precision. The scores written to the output will be output to this many decimal places. For example, a value of 2 will result in values
                                         output in the form '23.45', and a value of 3 will result in output of the form '23.449.
                                         
  -o [ --output ] arg (=norm_output.tsv) The name of the file to write output to. The output is formatted in the same way the input 'peptide_scores' are formatted, i.e. a score matrix with 
                                         samples on the columns and scores for a certain peptide on the rows. The score for each peptide in the output has been normalized in the manner 
                                         specified.
```
