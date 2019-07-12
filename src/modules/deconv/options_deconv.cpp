#include <algorithm>
#include <sstream>

#include "options_deconv.h"

std::string options_deconv::get_arguments()
{
    std::ostringstream str_stream;

    // returns the 'string' representation of a bool
    auto bool_str = []( bool opt ) -> std::string
        { return opt == true ? "true" : "false"; };

    str_stream << "--linked           " << linked_fname  << "\n" << 
                  " --threshold        " << threshold <<  "\n" << 
                  " --enriched         " << enriched_fname <<  "\n" << 
                  " --single_threaded  " <<  bool_str( single_threaded ) <<  "\n" << 
                  " --output           " << output_fname <<  "\n" << 
                  "\n";

    return str_stream.str();
}
