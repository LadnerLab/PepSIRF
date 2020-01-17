### Module p_enrich
```
PepSIRF: Peptide-based Serological Immune Response Framework species Paired (Duplicate) Enrichment module. This module determines which probes in samples done in duplicate are enriched, as determined by various numerical thresholds. For this module, thresholds are given as a comma-delimited pair. In order for the threshold to be met, each individual threshold in the pair must be met by at least one of the probes under consideration, independent of order. Note that a probe must meet each specified numeric threshold in order to be considered enriched.
:
  -h [ --help ]                     Produce help message and exit.
                                    
  -s [ --samples ] arg              The name of the file containing sample 
                                    pairs, denoting which samples are in 
                                    duplicate. This file must be tab-delimited 
                                    with one pair or samples per line.
                                    
  --zscores arg                     A matrix containing the zscores of each 
                                    probe in every sample. This should be in 
                                    the format output by the zscore module, 
                                    with probes on the rows and sample names on
                                    the columns.
                                    
  --zscore_constraint arg           Comma-separated values zscores should be 
                                    constrained by. These scores will be 
                                    evaluated as specified in the module 
                                    description. Example: '--zscore_constraint 
                                    3.5,4.0'. 
                                    
  --norm_scores arg                 A matrix containing normalized scores for 
                                    each probe in each sample.
                                    
  --norm_score_constraint arg       Comma-separated values normalized scores 
                                    should be constrained by. These scores will
                                    be evaluated as specified in the module 
                                    description. Example: '--norm_score_constra
                                    int 5.5,6.0'. 
                                    
  --raw_scores arg                  Optionally, a raw count matrix can be 
                                    included. This matrix must contain the raw 
                                    counts of each probe. If included, 
                                    'min_raw_count' must also be specified.
                                    
  --raw_score_constraint arg (=0,0) The minimum raw count a sample can have for
                                    all of its peptides in order for any of the
                                    probes in that sample to be considered 
                                    enriched. The sum of each probe's raw count
                                    in each sample must be at least either of 
                                    these values in order for the sample to be 
                                    considered.
                                    
  --outfile_suffix arg              Suffix to add to the names of the samples 
                                    written to output. For example, 
                                    '_enriched.txt' can be used. By default, no
                                    suffix is used.
                                    
  -j [ --join_on ] arg (=~)         A character or string to join output sample
                                    names on. For a sample pair of samples A 
                                    and B, the resulting file will have the 
                                    name 'A~B' if this flag is not given. 
                                    Otherwise, the given value will be used.
                                    
  -o [ --output ] arg (=paired)     Name of the directory to write output files
                                    to. Each sample with at least one enriched 
                                    peptide will receive a file in the output 
                                    directory.
```
