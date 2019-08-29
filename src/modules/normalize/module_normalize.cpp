#include <iostream>
#include <omp.h>
#include <algorithm>
#include <fstream>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/join.hpp>
#include <boost/range/adaptor/transformed.hpp>
#include <cstdlib>

#include "module_normalize.h"
#include "options_normalize.h"

module_normalize::module_normalize()
{
    name = "Normalize";
}

std::string module_normalize::get_name()
{
    return name;
}

void module_normalize::run( options *opts )
{
    options_normalize *n_opts = (options_normalize*) opts;

    std::string scores_fname = n_opts->peptide_scores_fname;

    omp_set_num_threads( n_opts->num_threads );

    peptide_score_data_sample_major original_scores;

    parse_peptide_scores( original_scores, scores_fname );

    write_peptide_scores( "pep_out.tsv", original_scores );

}

void module_normalize::parse_peptide_scores( peptide_score_data_sample_major& dest,
                                             std::string ifname
                                           )
{
    std::string line;
    std::ifstream in_file( ifname, std::ios_base::in );

    bool is_header = true; // first line is the header
    std::vector<std::string> split_line;

    while( std::getline( in_file, line ).good() )
        {
            boost::trim_right( line );
            boost::split( split_line, line, boost::is_any_of( "\t" ) );

            if( is_header )
                {
                    // copy all but the header (Sequence name) to the destination
                    std::for_each( split_line.begin() + 1, split_line.end(),
                                   [&]( const std::string &item )
                                   {  dest.sample_names.push_back( item ); }
                                 );

                    is_header = false;
                }
            else
                {
                    dest.pep_names.push_back( split_line[ 0 ] );
                    dest.scores.emplace_back( std::vector<std::size_t>() );

                    std::for_each( split_line.begin() + 1, split_line.end(),
                                   [&]( const std::string& item )
                                   { ( dest.scores.end() - 1)->push_back( std::atol( item.c_str() ) ); }
                                 );
                }
        }
}


void module_normalize::write_peptide_scores( std::string dest_fname,
                                             peptide_score_data_sample_major& data
                                           )
{
    std::ofstream out_file( dest_fname, std::ios_base::out );

    out_file << "Sequence name\t";

    out_file << boost::algorithm::join( data.sample_names, "\t" ) << "\n";

    for( std::size_t index = 0; index < data.pep_names.size(); ++index )
        {
            out_file << data.pep_names[ index ] << "\t" <<
                boost::algorithm::join( data.scores[ index ] |
                                        boost::adaptors::transformed(
                                                                     []( std::size_t s )
                                                                     { return std::to_string( s ); }
                                                                     ),
                                         "\t" ) << "\n";
        }
}
