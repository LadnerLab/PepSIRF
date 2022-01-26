#include "options_qc.h"
#include "options_demux.h"
#include <sstream>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/lexical_cast.hpp>


std::string options_qc::get_arguments()
{
    std::ostringstream str_stream;

    str_stream << "--index                    " << idx_fname << "\n" << 
                  " --index1                   " << tup_to_string( index1_data ) << "\n" <<
                  " --index2                   " << tup_to_string( index2_data ) << "\n" <<
                  " --samplelist               " << samplelist_fname << "\n"
                  "\n";
}