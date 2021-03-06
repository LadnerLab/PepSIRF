#### Module Deconv
```
PepSIRF: Peptide-based Serological Immune Response Framework species deconvolution module. This module has two modes: batch and singular. In batch mode, the input given to '--enriched' is a directory containing files of enriched peptides for each sample. Output is a directory where the output deconvolution reports will be written. In singular mode, both '--enriched' and '--output' are treated as files. The chosen mode is determined by the type of the argument to '--enriched'. If a directory is specified, batch mode will be used. If a file is specified, singular mode will be used.
:
  -h [ --help ]                         Produce help message
                                        
  -l [ --linked ] arg                   Name of file containing peptide to 
                                        species linkages. 
                                        
  -t [ --threshold ] arg                Threshold number of peptides for a 
                                        species to be considered. 
                                        
  --outfile_suffix arg                  Used for batch mode only. When 
                                        specified, each file written to the 
                                        output will have this suffix.
                                        
  -o [ --output ] arg (=deconv_output.tsv)
                                        Name of the directory or file to write 
                                        output to. Output will be in the form 
                                        of either one file or a directory 
                                        containing files, each a tab-delimited 
                                        file with a header. Each entry will be 
                                        of the form:
                                        species_id\tcount
                                         Note that in batch mode, the default 
                                        output directory will be 
                                        'deconv_output'.
                                        
  --scores_per_round arg                Name of directory to write 
                                        counts/scores to after every round. If 
                                        included, 
                                        the counts and scores for all remaining
                                        species will be written after every 
                                        round. 
                                        Filenames will be written in the format
                                        '$dir/round_x', where x is the round 
                                        number. 
                                        The original scores will be written to 
                                        '$dir/round_0'. A new file will be 
                                        written to the 
                                        directory after each subsequent round. 
                                        If this flag is included 
                                        and the specified directory exists, the
                                        program will exit with an error. 
                                        
  --single_threaded                     By default this module uses two 
                                        threads. Include this option with no 
                                        arguments if you only want  one thread 
                                        to be used. 
                                        
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
                                        shares a kmer with. 
                                        
  --summation_scoring                   Include this flag (without any 
                                        arguments) if you want summation 
                                        scoring to be used instead of 
                                        fractional or integer scoring. For 
                                        summation scoring, the --linked file 
                                        passed must be of the form created by 
                                        the link module. This means a file of 
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
                                        shares a kmer with. 
                                        
  -e [ --enriched ] arg                 Name of a directory containing files, 
                                        or a single file containing the names 
                                        of enriched peptides, one per line. 
                                        Each file in this file should have a 
                                        corresponding entry in the file 
                                        provided by the --linked option. 
                                        
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
                                        
  --peptide_assignment_map arg          If specified, a map detailing which 
                                        peptides were assigned to which species
                                        will be written. If deconv is used in 
                                        batch mode, this will be used as a 
                                        directory name for the peptide maps to 
                                        be stored. This map will be a 
                                        tab-delimited file with the  first 
                                        column peptide names, the second column
                                        is a comma-separated list of species 
                                        the peptide was assigned to. The third 
                                        column will be a list of the species 
                                        the peptide originally shared a kmer 
                                        with. 
                                        Note that the second column will only 
                                        contain multiple values in the event of
                                        a tie. 
                                        
  --mapfile_suffix arg                  In batch mode, add a suffix to the 
                                        filenames written to the 
                                        peptide_assignment_map directory. 
  --score_tie_threshold arg (=0)        Threshold for two species to be 
                                        evaluated as a tie. Note that this 
                                        value can be either an integer or a 
                                        ratio that is in (0,1). When provided 
                                        as an integer this value dictates the 
                                        difference in score that is allowed for
                                        two species to be considered. For 
                                        example, if this flag is provided with 
                                        the value 0, then two or more species 
                                        must have the exact same score to be 
                                        tied. If this flag is provided with the
                                        value 4, then the scores of species 
                                        must be no greater than 4 to be 
                                        considered tied. So if species 1has a 
                                        score of 5, and species has a score 
                                        anywhere between the integer values in 
                                        [1,9], then these species will be 
                                        considered tied, and their tie will be 
                                        evaluated as dicated by the tie 
                                        evaluation strategy provided.If the 
                                        argument provided to this flag is in 
                                        (0, 1), then a species must have at 
                                        least this percentage of the species 
                                        with the maximum score to be tied. So 
                                        if species 1 has the highest score with
                                        a score of 9, and species 2 has a score
                                        of 5, then this flag must be provided 
                                        with value >= 4/5 = 0.8 for the species
                                        to be considered tied. Note that any 
                                        values provided to this flag that are 
                                        in the set { x: x >= 1 } - Z, where Z 
                                        is the set of integers, will result in 
                                        an error. So 4.45 is not a valid value,
                                        but both 4 and 0.45 are.  
                                        
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
                                        id. 
                                        
  --score_overlap_threshold arg (=1)    Once two species have been found to be 
                                        within 'score_tie_threshold' number of 
                                        peptides of one another, they are then 
                                        evaluated as a tie. For a two-way tie 
                                        where integer tie evaluation is used, 
                                        if the species share more than 
                                        score_overlap_threshold number of 
                                        peptides, then they are both reported. 
                                        An example value is 10. For ratio tie 
                                        evaluation, which is used when this 
                                        argument is provided with a value in 
                                        the interval (0,1), two species must 
                                        share at leat this amount of peptides 
                                        with each other. For example, suppose 
                                        species 1 shares 0.5 of its peptides 
                                        with species 2, but species 2 only 
                                        shares 0.1 of its peptides with species
                                        1. To use integer tie evaluation, where
                                        species must share an integer number of
                                        peptides, not a ratio of their total 
                                        peptides, provide this argument with a 
                                        value in the interval [1, inf). These 
                                        two will only be reported together if 
                                        score_overlap_threshold <= 0.1.  
```
