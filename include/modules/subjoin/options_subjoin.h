#ifndef OPTIONS_SUBJOIN_HH_INCLUDED
#define OPTIONS_SUBJOIN_HH_INCLUDED
#include <vector>
#include <string>

#include "options.h"

/**
 * Options for the subjoin module of the 
 * PepSIRF package.
 **/
class options_subjoin : public options
{
 public:

    options_subjoin(); 
    /**
     * Returns a string of the arguments provided 
     * to the module.
     **/
    std::string get_arguments();

    std::vector<std::pair<std::string,std::string>> matrix_name_pairs;

    /**
     * The name of the file to write output to.
     **/
    std::string out_matrix_fname;

    /**
     * Boolean option to determine whether 
     * sample names or peptide names should be
     * used.
     **/
    bool use_sample_names;
};

#endif // OPTIONS_SUBJOIN_HH_INCLUDED
