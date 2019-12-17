#include <sstream>

#include "options_subjoin.h"

std::string options_subjoin::get_arguments()
{
    std::ostringstream str_stream;

    str_stream << "--names_list " << names_list_fname << "\n " <<
                  "--scores     " << in_matrix_fname << "\n " <<
                  "--output     " << out_matrix_fname << "\n " << 
                  "\n";

    return str_stream.str();
}
