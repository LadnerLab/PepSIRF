#include "options_bin.h"
#include <sstream>

options_bin::options_bin()
    : DEFAULT_MIN_BIN_SIZE{ 300 },
      DEFAULT_ROUNDING_FACTOR{ 0 } {}

std::string options_bin::get_arguments()
{
    std::ostringstream str_stream;

    str_stream << "--scores   " << input_scores_fname << "\n " 
               << "--bin_size " << min_bin_size << "\n "
               << "--round_to " << rounding_factor << "\n "
               << "--output   " << output_bins_fname << "\n"
        ;

    return str_stream.str();
};


