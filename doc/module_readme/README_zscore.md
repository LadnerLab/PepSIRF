#### Module Zscore
```
PepSIRF: Peptide-based Serological Immune Response Framework zscore module:
  -h [ --help ]                         Produce help message
                                        
  -s [ --scores ] arg                   Name of the file to use as input. 
                                        Should be a score matrix in the format 
                                        as output by the demux or subjoin 
                                        modules.
                                        
  -b [ --bins ] arg                     Name of the file containing bins, one 
                                        bin per line. Each bin contains a 
                                        tab-separated list of peptide names.
                                        
  --trim arg (=2.5)                     Percentile of lowest and highest counts
                                        within a bin to ignore when calculating
                                        the mean and standard deviation. This 
                                        value must be in the range 
                                        [0.00,100.0].
                                        
  -o [ --output ] arg (=zscore_output.tsv)
                                        Name of the file to write output to. In
                                        this file, each peptide will be written
                                        with its z-score within each sample.
                                        
  --nan_report arg                      Name of the file to write out 
                                        information regarding peptides that are
                                        given a zscore of 'nan'. This can 
                                        happen in one of two ways when all of 
                                        the peptides in a bin have a score of 
                                        zero. This will be a tab-delimited 
                                        file, with three columns per line. The 
                                        first column will contain the name of 
                                        the probe, the second will be the name 
                                        of the sample, and the third the bin 
                                        number of the probe. This bin number 
                                        corresponds to the line in the bins 
                                        file the probe was found in.
                                        
  -t [ --num_threads ] arg (=2)         The number of threads to use for 
                                        analyses.
```
