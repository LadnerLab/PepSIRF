#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/lexical_cast.hpp>
#include <sstream>
#include "options_demux.h"


std::string options_demux::get_arguments()
{
    std::ostringstream str_stream;

    str_stream << "--input_r1         " << input_r1_fname << "\n" <<
                  " --input_r2         " << input_r2_fname << "\n" << 
                  " --output           " << output_fname << "\n" <<
                  " --aa_counts        " << aggregate_fname << "\n" << 
                  " --index            " << index_fname << "\n" << 
                  " --samplelist       " << samplelist_fname << "\n"
                  " --library          " << library_fname << "\n" <<
                  " --read_per_loop    " << read_per_loop << "\n" <<
                  " --num_threads      " << num_threads << "\n" <<
                  " --seq              " << tup_to_string( seq_data ) << "\n" <<
                  " --f_index_data     " << tup_to_string( f_index_data ) << "\n" <<
                  " --r_index_data     " << tup_to_string( r_index_data ) << "\n" <<
                  " --phred_base       " << phred_base << "\n" <<
                  " --phred_min_score  " << min_phred_score << "\n" <<
                  " --concatemer       " << concatemer << "\n" << 
                  "\n";

    return str_stream.str();
}

std::string options_demux::tup_to_string( std::tuple<std::size_t, std::size_t, std::size_t>& data )
{
    std::ostringstream st;

    st << std::get<0>( data ) << "," << std::get<1>( data ) << "," << std::get<2>( data );
    return st.str();

}

void options_demux::set_info( std::tuple<std::size_t, std::size_t, std::size_t>
                              options_demux:: *member, std::string info
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
