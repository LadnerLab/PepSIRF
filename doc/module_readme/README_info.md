### Module Info
```
PepSIRF: Peptide-based Serological Immune Response Framework info module:
  -h [ --help ]         Produce help message and exit.
                        This module is used to gather information about a score
                        matrix. By default, the number of samples and peptides 
                        in the matrix will be output. Additional flags may be 
                        used to perform different analyses. For each input flag
                        included, one output file will be written.
                        
  -i [ --input ] arg    An input score matrix to gather information from.
                        
  --get_samples arg     Name of the file to write sample names to. Output will 
                        be in the form of a file with no header, one sample 
                        name per line.
                        
  --get_probes arg      Name of the file to write probe names to. Output will 
                        be in the form of a file with no header, one probe name
                        per line.
                        
  --col_sums arg        Name of the file to write the sum of column scores to. 
                        Output will be a tab-delimited file with a header. The 
                        first entry in each column will be the name of the 
                        sample, and the second will be the sum of the scores 
                        each peptide had for the sample.
```
