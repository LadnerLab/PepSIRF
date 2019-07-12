#include <algorithm>
#include <sstream>

#include "options_deconv.h"

std::string options_deconv::get_arguments()
{
    std::ostringstream str_stream;

    str_stream << "--linked    " << linked_fname  << "\n" << 
                  " --threshold " << threshold <<  "\n" << 
                  " --enriched  " << enriched_fname <<  "\n" << 
                  " --output    " << output_fname <<  "\n" << 
                  "\n";

    return str_stream.str();
}
