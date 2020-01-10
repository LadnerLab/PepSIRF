### Module s_enrich
```
PepSIRF: Peptide-based Serological Immune Response Framework species Single-Replicate Enrichment module:
  -h [ --help ]                 Produce help message and exit.
                                This module determines which probes in a sample
                                are enriched, as determined by various 
                                numerical thresholds that. Note that a probe 
                                must meet each specified numeric threshold in 
                                order to be considered enriched.
                                
  --zscores arg                 A matrix containing the zscores of each probe 
                                in every sample. This should be in the format 
                                output by the zscore module, with probes on the
                                rows and sample names on the columns.
                                
  --min_zscore arg              The minimum zscore a probe must have in order 
                                to be considered enriched.
                                
  --norm_scores arg             A matrix containing normalized scores for each 
                                probe in each sample.
                                
  --min_norm_score arg          The minimum normalized score a probe must have 
                                in a sample in order to be considered enriched.
                                
  --raw_scores arg              Optionally, a raw count matrix can be included.
                                This matrix must contain the raw counts of each
                                probe. If included, 'min_raw_count' must also 
                                be specified.
                                
  --min_raw_scores arg (=0)     The minimum raw count a sample can have for all
                                of its peptides in order for any of the probes 
                                in that sample to be considered enriched. The 
                                sum of each probe's raw count in a sample must 
                                be at least this value in order for the sample 
                                to be considered.
                                
  --outfile_suffix arg          Suffix to add to the names of the samples 
                                written to output. For example, '_enriched.txt'
                                can be used. By default, no suffix is used.
                                
  -o [ --output ] arg (=single) Name of the directory to write output files to.
                                Each sample with at least one enriched peptide 
                                will receive a file in the output directory.
```
