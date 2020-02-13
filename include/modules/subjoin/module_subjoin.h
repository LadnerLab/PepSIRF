#ifndef MODULE_SUBJOIN_HH_DEFINED
#define MODULE_SUBJOIN_HH_DEFINED
#include<fstream>

#include "module.h"
#include "options_subjoin.h"
#include "peptide_scoring.h"
#include "omp_opt.h"

/**
 * A module to join multiple score matrices into one.
 * Filters can be included to determine which peptide names will be output.
 **/
class module_subjoin : public module
{
 public:
    module_subjoin();

    using name_replacement_list = std::vector<std::pair<std::string,std::string>>;

    /**
     * Run the subjoin module.
     **/
    void run( options *opts );

    /**
     * Parse a list of names from the input stream.
     * There should be one name per line in the stream.
     * @param dest The location to store names, 
     *        this method will store items found in file at 
     *        the end of the vector.
     * @param file A stream to read names from.
     * @returns a list of the names that should be replaced in 
     *          the output. Each entry in the list will be 
     *          a pair of strings, with the first being the original name
     *          and the second the new value.
     **/
    name_replacement_list
        parse_namelist( std::vector<std::string>& dest,
                        std::istream& file
                      );


    /**
     * Full outer-Join matrices using a particular resolution strategy for
     * duplicates. 
     * @param first peptide_score_data_sample_major containing the first matrix to join.
     * @param second peptide_score_data_sample_major containing the second matrix to join.
     * @param resolution_strategy The resolution strategy to use for duplicate sample names.
     *        The resolution strategies are as follows: 
     *        - ignore: duplicates are ignored, it is not specified which score is taken.
     *        - combine: Both duplicates are kept, and their scores are combined.
     *        - include: Both duplicates are kept, and each samplename is given a suffix with 
     *                   the name of the file it's from.
     * @returns A matrix whose values are the result of the join.
     *          The exact values will depend upon the join resolution strategy used.
     **/
    labeled_matrix<double,std::string>
        join_with_resolve_strategy( peptide_score_data_sample_major first,
                                    peptide_score_data_sample_major second,
                                    evaluation_strategy::duplicate_resolution_strategy
                                    resolution_strategy
                                  ) const ;

    std::string get_name();

};

#endif // MODULE_SUBJOIN_HH_INCLUDED
