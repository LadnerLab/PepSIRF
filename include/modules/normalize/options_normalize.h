#ifndef OPTIONS_NORMALIZE_HH_INCLUDED 
#define OPTIONS_NORMALIZE_HH_INCLUDED 

#include <string>

#include "options.h"

class options_normalize : public options
{
 public:

    /**
     * Sets the name of this module to "norm"
     **/
    options_normalize();

    /**
     * Returns a string containing the arguments 
     * supplied to the program in a vertically-aligned
     * fashion. 
     **/
    std::string get_arguments();

    /**
     * Name of the file that contains the 
     * scores for peptides within a sample.
     **/
    std::string peptide_scores_fname;

    /**
     * The name of the file to write the module's
     * output to.
     **/
    std::string output_fname;

    bool col_sum_norm;
    bool size_factors_norm;


};

#endif // OPTIONS_NORMALIZE_HH_INCLUDED
