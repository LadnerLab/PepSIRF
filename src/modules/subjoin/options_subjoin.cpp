#include <sstream>
#include <algorithm>

#include "options_subjoin.h"

std::string options_subjoin::get_arguments()
{
    std::ostringstream str_stream;

    std::for_each( matrix_name_pairs.begin(),
                   matrix_name_pairs.end(),
                   [&]( const std::pair<std::string,std::string>& str )
                   {
                       str_stream << "--filter_scores "
                                  << str.first << "," << str.second
                                  << "\n ";
                   }
                 );
    str_stream << "--output        " << out_matrix_fname << "\n " << 
                  "\n";

    return str_stream.str();
}
