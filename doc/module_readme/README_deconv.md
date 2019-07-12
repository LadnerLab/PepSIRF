#### Module Deconv
```
PepSIRF: Peptide-based Serological Immune Response Framework species deconvolution module. 
:
  -h [ --help ]          Produce help message
                         
  -l [ --linked ] arg    Name of file containing peptide to species linkages.
                         
  -t [ --threshold ] arg Threshold number of peptides for a species to be 
                         considered.
                         
  -o [ --output ] arg    Name of the file to write output to. Output will be in
                         the form of a tab-delimited file with a header. Each 
                         entry will be of the form:
                         species_id\tcount
                         
  --single_threaded      By default this module uses two threads. Include this 
                         option with no arguments if you only want  one thread 
                         to be used.
                         
  --fractional_scoring   Use fractional instead of integer scoring. For integer
                         scoring the score of each species is defined by the 
                         number of peptides that share a 7mer with that 
                         species. For fractional scoring the score of each 
                         species is defined by 1/n for each peptide, where n is
                         the number of species a peptide shares a 7mer with. In
                         this method of scoring peptides with fewer species are
                         worth more.
                         
  -e [ --enriched ] arg  File containing the names of enriched peptides, one 
                         per line. Each file in this file should have a 
                         corresponding entry in the file provided by the 
                         --linked option.
```                        
