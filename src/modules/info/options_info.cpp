#include "options_info.h"
#include <sstream>

options_info::options_info() = default;
options_info::~options_info() = default;

std::string options_info::get_arguments()
{
    std::ostringstream stream;

    stream
        << "--input       " << in_fname << "\n " 
        << "--get_samples " << out_samples_fname << "\n " 
        << "--get_names   " << out_pep_names_fname << "\n " 
        << "--col_sums    " << out_col_sums_fname << "\n " 
        << "\n"
        ;

    return stream.str();
}
