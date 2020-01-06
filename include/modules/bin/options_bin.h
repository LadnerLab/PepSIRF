#ifndef OPTIONS_BIN_HH_INCLUDED
#define OPTIONS_BIN_HH_INCLUDED
#include <string>

#include "options.h"

class options_bin : public options
{
 public:
    options_bin();
    std::string get_arguments();

    /**
     * Input filename containing scores 
     * of peptides in each sample. 
     **/
    std::string input_scores_fname;

    /**
     * Name of the file to write output bins to
     **/
    std::string output_bins_fname;

    /**
     * The fewest peptides a bin can have.
     **/
    std::size_t min_bin_size;

    /**
     * Round to the nearest 1/ 10^(rounding_factor)
     **/
    std::size_t rounding_factor;

    const std::size_t DEFAULT_MIN_BIN_SIZE;
    const std::size_t DEFAULT_ROUNDING_FACTOR;

};

#endif // OPTIONS_BIN_HH_INCLUDED
