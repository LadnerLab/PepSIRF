#include <iostream>
#include <cstdio>
#include "omp_opt.h"
#include <algorithm>
#include <fstream>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/join.hpp>
#include <boost/range/adaptor/transformed.hpp>
#include <cstdlib>
#include <functional>
#include <numeric>

#include "module_normalize.h"
#include "options_normalize.h"
#include "time_keep.h"
#define ONE_MILLION 1000000

module_normalize::module_normalize()
{
    name = "Norm";
}

std::string module_normalize::get_name()
{
    return name;
}

void module_normalize::run( options *opts )
{
    options_normalize *n_opts = (options_normalize*) opts;
    time_keep::timer timer;

    timer.start();

    std::string scores_fname = n_opts->peptide_scores_fname;

    omp_set_num_threads( n_opts->num_threads );

    peptide_score_data_sample_major original_scores;

    parse_peptide_scores( original_scores, scores_fname );

    std::vector<double> norm_factors( original_scores.sample_names.size(), 0 );

    if( n_opts->size_factors_norm )
        {
            compute_size_factors( norm_factors, original_scores.scores );
        }
    else
        {
            get_sum( norm_factors, original_scores.scores );
            constant_factor_normalization( norm_factors, ONE_MILLION );
        }

    // normalize the counts
    normalize_counts( original_scores.scores, norm_factors );

    write_peptide_scores( n_opts->output_fname, original_scores );

    timer.stop();
    
    std::cout << "Took " << timer.get_elapsed() << " seconds.\n";

}

void module_normalize::parse_peptide_scores( peptide_score_data_sample_major& dest,
                                             std::string ifname
                                           )
{
    std::string line;
    std::ifstream in_file( ifname, std::ios_base::in );

    std::vector<std::string> lines_from_file;
    std::vector<std::string> split_line;

    while( std::getline( in_file, line ).good() )
        {
            boost::trim_right( line );
            lines_from_file.emplace_back( line );
        }

    std::size_t line_count = lines_from_file.size();
    std::size_t num_peptides = line_count - 1; // -1 because of the header
    std::size_t sample_count = std::count( lines_from_file[ 0 ].begin(),
                                           lines_from_file[ 0 ].end(),
                                           '\t'
                                         );
    dest.scores.reserve( num_peptides );
    dest.pep_names.reserve( num_peptides );
    dest.sample_names.reserve( sample_count );

    std::size_t index = 0;

    for( index = 0; index < num_peptides; ++index )
        {
            dest.pep_names.emplace_back( std::string( "" ) );
        }

    for( index = 0; index < num_peptides; ++index )
        {
            dest.scores.emplace_back( std::vector<double>( sample_count, 0 ) );
            dest.scores.back().reserve( sample_count );
        }

    // save the sample names to the samplenames vector
    boost::split( split_line, lines_from_file[ 0 ], boost::is_any_of( "\t" ) );
    std::for_each( split_line.begin() + 1, split_line.end(),
                   [&]( const std::string &item )
                   {
                       dest.sample_names.push_back( item );
                   }
                 );

    for( index = 1; index < lines_from_file.size(); ++index )
        {
            std::size_t assign_index = index - 1;
            // get the name of the peptide
            const std::string& my_str = lines_from_file[ index ];
            std::size_t pos = my_str.find( '\t' );

            dest.pep_names[ assign_index ] = my_str.substr( 0, pos );

            std::size_t begin_pos = pos + 1;
            
            // for j = 0 -> num_samples grab the number
            for( std::size_t count_index = 0; count_index < sample_count; ++count_index )
                {
                    std::size_t end_pos   = my_str.find( '\t', begin_pos );
                    double val = std::strtod( my_str.substr( begin_pos, end_pos - begin_pos ).c_str(),
                                     nullptr
                                   );
                    dest.scores[ assign_index ][ count_index ] = val;
                    begin_pos = end_pos + 1;
                        
                }
        }

}


void module_normalize::write_peptide_scores( std::string dest_fname,
                                             peptide_score_data_sample_major& data
                                           )
{
    std::ofstream out_file( dest_fname, std::ios_base::out );
    char digits[ 30 ];

    out_file << "Sequence name\t";

    out_file << boost::algorithm::join( data.sample_names, "\t" ) << "\n";

    for( std::size_t index = 0; index < data.pep_names.size(); ++index )
        {
            out_file << data.pep_names[ index ] << "\t";
                        
            std::size_t inner_index = 0;

            for( inner_index = 0; inner_index < data.sample_names.size(); ++inner_index )
                {
                    std::sprintf( digits, "%.2f", data.scores[ index ][ inner_index ] );
                    out_file << digits;

                    if( inner_index < data.sample_names.size() - 1 )
                        {
                            out_file << "\t";
                        }
                }
            out_file << "\n";
        }
}

void module_normalize::constant_factor_normalization( std::vector<double>& cols,
                                                      const std::size_t factor
                                                    )
{
    for( std::size_t index = 0; index < cols.size(); ++index )
        {
            cols[ index ] /= factor;
        }
}

void module_normalize::get_sum( std::vector<double>& dest,
                                std::vector<std::vector<double>>& src
                                )
{
    std::size_t index = 0;

    for( index = 0; index < src.size(); ++index )
        {
            std::size_t inner_index = 0;
            for( inner_index = 0; inner_index < src[ index ].size(); ++inner_index )
                {
                    dest[ inner_index ] += src[ index ][ inner_index ];
                }
            
        }
}

void module_normalize::normalize_counts( std::vector<std::vector<double>>&
                                         original_scores,
                                         const std::vector<double>& norm_factors
                                       )
{
    for( std::size_t index = 0; index < original_scores.size(); ++index )
        {
            for( std::size_t inner_index = 0; inner_index < original_scores[ index ].size(); ++inner_index )
                {
                    original_scores[ index ][ inner_index ] /= norm_factors[ inner_index ];
                }
        }
}

void module_normalize::compute_size_factors( std::vector<double>& size_factors,
                                             const std::vector<std::vector<double>>& data
                                           )
{
 //    auto median = []( std::vector<double>& v ) -> double
 //        {
 // size_t n = v.size() / 2;
 //  std::nth_element(v.begin(), v.begin()+n, v.end());
 //  double vn = v[n];
 //  if(v.size()%2 == 1)
 //  {
 //    return vn;
 //  }else
 //  {
 //    std::nth_element(v.begin(), v.begin()+n-1, v.end());
 //    return 0.5*(vn+v[n-1]);
 //  }
 //        };

    std::vector<double> row;
    std::vector<std::vector<double>> geom_mean_factors;
    geom_mean_factors.reserve( data[ 0 ].size() );

    std::size_t index = 0;

    for( index = 0; index < data[ 0 ].size(); ++index )
        {
            geom_mean_factors.emplace_back( std::vector<double>( data.size(), 0 ) ) ;
        }

    for( index = 0; index < data.size(); ++index )
        {
            double g_mean = geom_mean( data[ index ] );

            for( std::size_t inner_index = 0; inner_index < data[ index ].size(); ++inner_index )
                {
                    // transpose the data here so computing the medians is easier
                    geom_mean_factors[ inner_index ][ index ] = data[ index ][ inner_index ] / g_mean;
                }
        }

    for( index = 0; index < geom_mean_factors.size(); ++index )
        {
            std::sort( geom_mean_factors[ index ].begin(),
                       geom_mean_factors[ index ].end()
                     );

            auto begin = std::upper_bound( geom_mean_factors[ index ].end(),
                                           geom_mean_factors[ index ].end(),
                                           0
                                         );
            int median_loc = (geom_mean_factors[ index ].size()
                              - std::distance( begin, geom_mean_factors[ index ].end() )) / 2;
            size_factors.push_back( geom_mean_factors[ index ][ median_loc ] );
        }

}

