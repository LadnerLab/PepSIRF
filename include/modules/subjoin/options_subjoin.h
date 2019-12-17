#ifndef OPTIONS_SUBJOIN_HH_INCLUDED
#define OPTIONS_SUBJOIN_HH_INCLUDED
#include "options.h"

/**
 * Options for the subjoin module of the 
 * PepSIRF package.
 **/
class options_subjoin : public options
{
 public:

    /**
     * Returns a string of the arguments provided 
     * to the module.
     **/
    std::string get_arguments();

    /**
     * The filename of the file containing 
     * the names of peptides to keep.
     **/
    std::string names_list_fname;

    /**
     * The name of the matrix file to be 
     * filtered.
     **/
    std::string in_matrix_fname;

    /**
     * The name of the file to write output to.
     **/
    std::string out_matrix_fname;
};

#endif // OPTIONS_SUBJOIN_HH_INCLUDED
