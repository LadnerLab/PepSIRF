#include <iostream>
#include <cstdio>
#include "omp_opt.h"
#include <algorithm>
#include <fstream>
#include <cstdlib>
#include <functional>
#include <numeric>
#include <iomanip>
#include <boost/algorithm/string.hpp>

#include "module_normalize.h"
#include "options_normalize.h"
#include "time_keep.h"
#include "stats.h"
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
    peptide_score_data_sample_major neg_scores;
    peptide_scoring::parse_peptide_scores( original_scores, scores_fname );
    original_scores.scores = original_scores.scores.transpose();
    std::vector<std::string> neg_filter;

    if( !(n_opts->neg_control).empty() )
        {
            peptide_scoring::parse_peptide_scores( neg_scores, n_opts->neg_control );

            if( !(n_opts->neg_names).empty() )
                {
                    boost::split( neg_filter, neg_scores.sample_names, boost::is_any_of( "," ) );
                }
            else if( !(n_opts->neg_id).empty() )
                {
                    for( const auto& sample : neg_scores.sample_names )
                        {
                            if( sample.rfind( n_opts->neg_id, 0 ) != std::string::npos )
                                {
                                    neg_filter.emplace_back( sample );
                                }
                        }
                }



        }

    std::size_t sample_size;

    if( !(n_opts->neg_control).empty() )
        {
            sample_size = neg_scores.sample_names.size();
        }
    else
        {
            sample_size = original_scores.sample_names.size();
        }

    std::vector<double> norm_factors( sample_size, 0 );

    if( n_opts->approach == "size_factors" )
        {
            compute_size_factors( norm_factors, original_scores.scores );
        }
    else if( n_opts->approach == "neg_diff" )
        {
            // get ratio, pass in norm_factors, the original scores, the neg scores,
            // and the mean average of the negative scores.
        }
    else if( n_opts->approach == "diff_ratio" )
        {

        }
    else // Column sum default
        {
            get_sum( norm_factors, original_scores.scores );
            constant_factor_normalization( norm_factors, ONE_MILLION );
        }

    // normalize the counts
    normalize_counts( original_scores.scores, norm_factors );

    original_scores.scores = original_scores.scores.transpose();

    std::ofstream output_file( n_opts->output_fname, std::ios_base::out );
    output_file
        << std::fixed
        << std::showpoint
        << std::setprecision( n_opts->precision_digits );

    peptide_scoring::write_peptide_scores( output_file, original_scores );

    timer.stop();

    std::cout << "Took " << timer.get_elapsed() << " seconds.\n";

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

void module_normalize::get_ratio( std::vector< )

void module_normalize::get_sum( std::vector<double>& dest,
                                matrix<double>& src
                              )
{
    std::size_t index = 0;

    for( index = 0; index < src.nrows(); ++index )
        {
            std::size_t inner_index = 0;
            for( inner_index = 0; inner_index < src.ncols(); ++inner_index )
                {
                    dest[ inner_index ] += src( index, inner_index );
                }

        }
}

void module_normalize::normalize_counts( matrix<double>&
                                         original_scores,
                                         const std::vector<double>& norm_factors
                                       )
{
    for( std::size_t index = 0; index < original_scores.nrows(); ++index )
        {
            for( std::size_t inner_index = 0; inner_index < original_scores.ncols(); ++inner_index )
                {
                    original_scores( index, inner_index ) /= norm_factors[ inner_index ];
                }
        }
}

void module_normalize::compute_size_factors( std::vector<double>& size_factors,
                                             const matrix<double>& data
                                           )
{
    std::vector<double> row;
    std::vector<std::vector<double>> geom_mean_factors;
    geom_mean_factors.reserve( data.ncols() );

    std::size_t index = 0;

    for( index = 0; index < data.nrows(); ++index )
        {
            geom_mean_factors.emplace_back( std::vector<double>( data.nrows(), 0 ) ) ;
        }

    for( index = 0; index < data.nrows(); ++index )
        {
            double g_mean = stats::geom_mean( data.row_begin( index ),
                                              data.row_end( index )
                                            );

            for( std::size_t inner_index = 0; inner_index < data.ncols(); ++inner_index )
                {
                    // transpose the data here so computing the medians is easier
                    geom_mean_factors[ inner_index ][ index ] = data( index, inner_index ) / g_mean;
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

