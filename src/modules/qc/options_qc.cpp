#include "options_qc.h"
#include <sstream>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/lexical_cast.hpp>


std::string options_qc::get_arguments()
{
    std::ostringstream str_stream;

    str_stream << " --fif                    " << idx_fname << "\n" << 
                  " --input_r1                   " << input_r1_fname << "\n" <<
                  " --input_r2                   " << input_r2_fname << "\n" <<
                  " --samplelist               " << samplelist_fname << "\n" <<
                  " --outfile                   " << output_fname << "\n" <<
                  " --sname                     " << samplename << "\n" <<
                  " --index1                    " << tup_to_string( index1_data ) << "\n" <<
                  " --index2                    " << tup_to_string( index2_data ) << "\n" <<
                  "\n";

    return str_stream.str();
}

std::string options_qc::tup_to_string( std::tuple<std::size_t, std::size_t, std::size_t>& data )
{
    std::ostringstream st;

    st << std::get<0>( data ) << "," << std::get<1>( data ) << "," << std::get<2>( data );
    return st.str();
}

void options_qc::set_info( std::tuple<std::size_t, std::size_t, std::size_t>
                              options_qc:: *member, std::string info
                            )
{
    const char *COMMA             = ",";
    const int  NUM_REQUIRED_ITEMS = 3;

    std::vector<std::string> matches;
    boost::split( matches, info, boost::is_any_of( COMMA ) );

    auto cast_fn = [=]( std::string m ) -> std::size_t
        { return boost::lexical_cast<std::size_t>( m ); };

    if( matches.size() != NUM_REQUIRED_ITEMS )
        {
            throw std::runtime_error( "Incorrect number of comma-separated values "
                                      "provided. Please try ./pepsirf demux -h for help"
                                    );
        }

    auto tup = std::make_tuple( cast_fn( matches[ 0 ] ),
                                cast_fn( matches[ 1 ] ),
                                cast_fn( matches[ 2 ] )
                              );
    this->*member = tup;
}