#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/lexical_cast.hpp>
#include <sstream>
#include <vector>
#include "options_demux.h"
#include "logger.h"


std::string options_demux::get_arguments()
{
    std::ostringstream str_stream;
    str_stream
        << "--input_r1              " << input_r1_fname << "\n"
        << "--input_r2              " << input_r2_fname << "\n"
        << "--output                " << output_fname << "\n"
        << "--aa_counts             " << aggregate_fname << "\n"
        << "--translate_aggregates  " << std::boolalpha << translation_aggregation << "\n"
        << "--index                 " << index_fname << "\n"
        << "--fif                   " << flexible_idx_fname << "\n"
        << "--samplelist            " << samplelist_fname << "\n"
        << "--sname                 " << samplename << "\n"
        << "--sindex                " << indexes << "\n"
        << "--library               " << library_fname << "\n"
        << "--read_per_loop         " << read_per_loop << "\n"
    	<< "--num_threads           " << num_threads << "\n"
        << "--include_toggle        " << pos_toggle << "\n"
        << "--seq                   " << tup_to_string( seq_data ) << "\n"
        << "--index1                " << tup_to_string( index1_data ) << "\n"
        << "--index2                " << tup_to_string( index2_data ) << "\n"
        << "--phred_base            " << phred_base << "\n"
        << "--phred_min_score       " << min_phred_score << "\n"
        << "--concatemer            " << concatemer << "\n"
        << "--diagnostic_info       " << diagnostic_fname << "\n"
        << "--fastq_output          " << fastq_out << "\n"
        << "--logfile               " << logfile << "\n"
        << "--replicate_info        " << replicate_info_fname << "\n"
        << "--unmapped_reads_output " << unmapped_reads_fname <<  "\n"
        << "--trunc_info_output     " << trunc_info_outdir <<  "\n\n";

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
            Log::error(
                "Incorrect number of comma-separated values provided."
                " Please try ./pepsirf demux -h for help"
            );
        }

    auto tup = std::make_tuple( cast_fn( matches[ 0 ] ),
                                cast_fn( matches[ 1 ] ),
                                cast_fn( matches[ 2 ] )
                              );
    this->*member = tup;
}
