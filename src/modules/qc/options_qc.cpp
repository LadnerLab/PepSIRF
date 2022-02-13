#include "options_qc.h"
#include <sstream>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/lexical_cast.hpp>


std::string options_qc::get_arguments()
{
    std::ostringstream str_stream;

    str_stream << " --index                    " << idx_fname << "\n" << 
                  " --input_r1                   " << tup_to_string( index1_data ) << "\n" <<
                  " --input_r2                   " << tup_to_string( index2_data ) << "\n" <<
                  " --samplelist               " << samplelist_fname << "\n" <<
                  " --sname                     " << samplename <<
                  " --outfile                   " << output_fname <<
                  "\n";

    return str_stream.str();
}

std::string options_qc::tup_to_string( std::tuple<std::size_t, std::size_t, std::size_t>& data )
{
    std::ostringstream st;

    st << std::get<0>( data ) << "," << std::get<1>( data ) << "," << std::get<2>( data );
    return st.str();
}