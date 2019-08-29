#include <iostream>
#include <omp.h>
#include <algorithm>
#include <fstream>
#include <boost/algorithm/string.hpp>
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

