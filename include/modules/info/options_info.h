#ifndef OPTIONS_INFO_HH_INCLUDED
#define OPTIONS_INFO_HH_INCLUDED

#include <string>
#include "options.h"

class options_info : public options
{

public:

    options_info();
    ~options_info();

    std::string get_arguments();

    /**
     * The name of the input score matrix.
     **/
    std::string in_fname;

    /**
     * The name of the file to write output 
     * sample names to.
     **/
    std::string out_samples_fname;

    /**
     * The name of the file to write output 
     * peptide names to.
     **/
    std::string out_pep_names_fname;

    /**
     * The name of the file to write 
     * output column sums to
     **/
    std::string out_col_sums_fname;
    
    /**
    * The name of the input replicate names file
    **/
    std::string in_replicates_fname;
    
    /**
     * The name of the input replicate names file
     **/
    std::string out_avgs_fname;
};

#endif // OPTIONS_INFO_HH_INCLUDED
