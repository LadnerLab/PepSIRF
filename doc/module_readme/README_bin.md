#### Module Deconv
```
PepSIRF: Peptide-based Serological Immune Response Framework bin module.
:
  -h [ --help ]                         Produce help message and exit.
                                        Typically, the subjoin module will be 
                                        used to specify samples that are 
                                        negative controls. 
                                        This module is usually used to create 
                                        bins from these, negative controls, 
                                        allowing the z-scores for a peptide 
                                        to be calculated from the scores of the
                                        peptides it shares a bin with.
                                        
  -s [ --scores ] arg                   The input file to parse scores from. 
                                        Bins will be created from the peptides 
                                        in these bins and the scores the 
                                        peptides have. Peptides with similar 
                                        scores will be binned together.
                                        
  -b [ --bin_size ] arg (=300)          The minimum number of peptides that can
                                        be placed in a bin. If a bin would be 
                                        created with fewer than bin_size 
                                        peptides, it will be combined with the 
                                        next bin until at least bin_size 
                                        peptides are found.
                                        
  -r [ --round_to ] arg (=0)            The 'rounding factor' for the scores 
                                        parsed from the score matrix. Scores 
                                        found in the matrix will be rounded to 
                                        the nearest 1/10^x for a rounding 
                                        factor x. For example, a rounding 
                                        factor of 0 will result in rounding to 
                                        the nearest integer, while a rounding 
                                        factor of 1 will result in rounding to 
                                        the nearest tenth.
                                        
  -o [ --output ] arg (=bin_output.tsv) Name of the file to write the bins to, 
                                        one per line. Each bin will be a 
                                        tab-delimited list of the names of the 
                                        peptides in the bin.
```                                        
