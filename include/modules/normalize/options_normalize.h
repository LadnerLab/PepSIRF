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

    std::string approach;
    /* remove these two variables -> approach changing
    to take a character to specify the normalization approach

    bool col_sum_norm;
    bool size_factors_norm;
    */
    std::string neg_control;

    std::string neg_names;

    std::string neg_id;

    /**
     * Number of digits of precision for the output.
     **/
    int precision_digits;


};

#endif // OPTIONS_NORMALIZE_HH_INCLUDED
