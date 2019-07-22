#### Module Deconv
```
PepSIRF: Peptide-based Serological Immune Response Framework species deconvolution module. 
This module has two different modes: scoring species and creating a linkage file. Each of these modes has its own set of arguments and parameters. The description of each module is followed by the mode for which this command is, in brackets. For example, if the description is followed by [scoring_species], then this argument is for the scoring species mode. Similarly, [create_linkage] is followed by linkage creation arguments. Arguments pertinent to both modes are followed by [scoring_species,create_linkage].
:
  -h [ --help ]                         Produce help message
                                        
  -l [ --linked ] arg                   Name of file containing peptide to 
                                        species linkages. [scoring_species]
                                        
  -t [ --threshold ] arg                Threshold number of peptides for a 
                                        species to be considered. 
                                        [scoring_species]
                                        
  -o [ --output ] arg (=deconv_output.tsv)
                                        Name of the file to write output to. 
                                        Output will be in the form of a 
                                        tab-delimited file with a header. Each 
                                        entry will be of the form:
                                        species_id\tcount
                                         [create_linkage,scoring_species]
                                        
  --single_threaded                     By default this module uses two 
                                        threads. Include this option with no 
                                        arguments if you only want  one thread 
                                        to be used. [create_linkage,scoring_spe
                                        cies]
                                        
  --fractional_scoring                  Use fractional instead of integer 
                                        scoring. For integer scoring the score 
                                        of each species is defined by the 
                                        number of peptides that share a 7mer 
                                        with that species. For fractional 
                                        scoring the score of each species is 
                                        defined by 1/n for each peptide, where 
                                        n is the number of species a peptide 
                                        shares a 7mer with. In this method of 
                                        scoring peptides with fewer species are
                                        worth more. Note that if neither this 
                                        flag nor --summation_scoring are 
                                        included, integer scoring will be used.
                                        In integer scoring each species is 
                                        scored by the number of peptides it 
                                        shares a kmer with. [scoring_species]
                                        
  --summation_scoring                   Include this flag (without any 
                                        arguments) if you want summation 
                                        scoring to be used instead of 
                                        fractional or integer scoring. For 
                                        summation scoring, the --linked file 
                                        passed must be of the form created by 
                                        --create_linkage. This means a file of 
                                        tab-delimited values, one per line. 
                                        Each line is of the form peptide_name 
                                        TAB id:score,id:score, and so on. 
                                        Undefined behavior will result if input
                                        is not in this format. For summation 
                                        scoring, each species is scored based 
                                        on the number of kmers it shares with 
                                        each peptide with which it shares a 
                                        kmer.
                                         For example, assume a line in the 
                                        --linked file looks like the following:
                                         
                                        peptide_1 TAB 123:4,543:8
                                        Both species '123' and '543' will 
                                        receive a score of 4 and 8 
                                        respectively.Note that if neither this 
                                        flag nor --summation_scoring are 
                                        included, integer scoring will be used.
                                        In integer scoring each species is 
                                        scored by the number of peptides it 
                                        shares a kmer with. [scoring_species]
                                        
  -e [ --enriched ] arg                 File containing the names of enriched 
                                        peptides, one per line. Each file in 
                                        this file should have a corresponding 
                                        entry in the file provided by the 
                                        --linked option. [scoring_species]
                                        
  --score_filtering                     Include this flag if you want filtering
                                        to be done by the score of each 
                                        species. Note that score is determined 
                                        by the different flags specifying how a
                                        species should be scored. This means 
                                        that any species whose score falls 
                                        below --threshold will be removed from 
                                        consideration. Note that for integer 
                                        scoring, both score filtering and count
                                        filtering are the same. If this flag is
                                        not included, then any species whose 
                                        count falls below --threshold will be 
                                        removed from consideration. Score 
                                        filtering is best suited for the 
                                        summation scoring algorithm. 
                                        [scoring_species]
                                        
  --create_linkage                      Boolean switch to create the linkage 
                                        file that is used as input for the 
                                        species deconvolution process. If this 
                                        option is included then 'protein_file' 
                                        and 'peptide_file' must also be 
                                        included. [create_linkage]
                                        
  --protein_file arg                    Name of fasta file containing protein 
                                        sequences from which a design was 
                                        created. [create_linkage]
                                        
  --peptide_file arg                    Name of fasta file containing aa 
                                        peptides that have been designed as 
                                        part of a library. [create_linkage]
                                        
  -k [ --kmer_size ] arg                Kmer size to use when creating the 
                                        linkage map. Entries in the linkage 
                                        file will contain peptides and the 
                                        species ids of species that share a 
                                        kmer with this peptide. For example, if
                                        k is 7 and there exists a line in the 
                                        linkage file of the form: 
                                         'peptide_1 TAB 455:12,423:10'
                                         then peptide_1 shares 12 7-mers with 
                                        the species with id '455', and 10 
                                        7-mers with the species that has id 
                                        423. [create_linkage]
                                        
  --id_name_map arg                     File containing mappings from taxonomic
                                        id to name. This file should be 
                                        formatted like the file 
                                        'rankedlineage.dmp' from NCBI. It is 
                                        recommended to either use this file or 
                                        a subset of this file that at least 
                                        contains the species ids of the 
                                        designed peptides. If included, the 
                                        output will contain a column denoting 
                                        the name of the species as well as the 
                                        id. [create_linkage]
                                        
```                        

