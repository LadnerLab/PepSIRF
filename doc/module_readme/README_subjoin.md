#### Module Subjoin
```
PepSIRF: Peptide-based Serological Immune Response Framework subjoin module:
  -h [ --help ]                         Produce help message
                                        
  --filter_scores arg                   Either comma-separated filenames (For 
                                        example: score_matrix.tsv,sample_names.
                                        txt ) or the name of a tab-delimited 
                                        file containing score_matrix and sample
                                        name list filename pairs, one per line.
                                        Each of these pairs must be a score 
                                        matrix and a file containing the names 
                                        of samples (or peptides, if specified) 
                                        to keep in the score matrix. The score 
                                        matrix should be of the format output 
                                        by the demux module, with sample names 
                                        on the columns and peptide names on the
                                        rows. The namelist must have one name 
                                        per line. To use multiple name lists 
                                        with multiple score matrices, include 
                                        this argument multiple times.
                                        
  --filter_peptide_names                Flag to include if the files input to 
                                        the filter_scores options should be 
                                        treated as peptide names instead of 
                                        sample names. With the inclusion of 
                                        this flag, the input files will be 
                                        filtered on peptide names (rows) 
                                        instead of sample names (column).
                                        
  --duplicate_evaluation arg (=include) Defines what should be done when sample
                                        or peptide names are not unique across 
                                        files being joined. Currently, three 
                                        different duplicate evaluation 
                                        strategies are used: 
                                         - combine: Combine (with addition) the
                                        values associated with peptide/sample 
                                        names. 
                                        
                                         - include: Include each duplicate, 
                                        adding a suffix to the duplicate 
                                        samplename detailing the file the 
                                        sample came from. 
                                        
                                         - ignore: Ignore the possibility of 
                                        duplicates. Behavior is undefined when 
                                        duplicates are  encountered in this 
                                        mode.
                                        
                                         Possible options include combine, 
                                        include, and ignore.
                                        
  --output arg (=subjoin_output.tsv)    The name of the file to write output 
                                        scores to. The output will be in the 
                                        form of the input, but with only the 
                                        specified values (samplenames or 
                                        peptides) found in the namelists. 
```
