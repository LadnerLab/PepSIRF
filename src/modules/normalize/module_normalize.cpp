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
    std::size_t sample_size;

    // Use negative control if provided
    if( !n_opts->neg_control.empty() )
        {
            peptide_scoring::parse_peptide_scores( neg_scores, n_opts->neg_control );
// This can be made into a function where the data sample major is passed (remove repetive code)
// Function : filter_neg_control_start_names
            // Filter negative scores by: string identifiable prefix or list of sample names
            if( !n_opts->neg_names.empty() && !n_opts->neg_id.empty() )
                {
                    if( !neg_scores.file_name.empty() )
                        {
                            for( const auto& sample : neg_scores.sample_names )
                                {
                                    if( sample.rfind( n_opts->neg_id, n_opts->neg_id.size() ) != std::string::npos
                                        || boost::contains( n_opts->neg_names, sample ) )
                                        {
                                            neg_filter.emplace_back( sample );
                                        }
                                }
                        }
                    else
                        {
                            for( const auto& sample : original_scores.sample_names )
                                {
                                    if( sample.rfind( n_opts->neg_id, n_opts->neg_id.size() ) != std::string::npos
                                        || boost::contains( n_opts->neg_names, sample ) )
                                        {
                                            neg_filter.emplace_back( sample );
                                        }
                                }
                        }
                }
            else if( !n_opts->neg_names.empty() )
                {
                    boost::split( neg_filter, n_opts->neg_names, boost::is_any_of( "," ) );
                }
            else if( !n_opts->neg_id.empty() )
                {
// Function : filter_neg_control_start
                    if( !neg_scores.file_name.empty() )
                        {
                            for( const auto& sample : neg_scores.sample_names )
                                {
                                    if( sample.rfind( n_opts->neg_id, n_opts->neg_id.size() ) != std::string::npos )
                                        {
                                            neg_filter.emplace_back( sample );
                                        }
                                }
                        }
                    else // if data matrix of sb/neg samples is not provided, use input matrix
                        {
                            for( const auto& sample : original_scores.sample_names )
                                {
                                    if( sample.rfind( n_opts->neg_id, n_opts->neg_id.size() ) != std::string::npos )
                                        {
                                            neg_filter.emplace_back( sample );
                                        }
                                }
                        }

                }
            else
                {
                    throw std::runtime_error( "Error: Must use approach for identifying negative controls. "
                                         "Either a negative id [--negative_id,-s] or [--negative_names,-n]." );
                }
        }

    if( !n_opts->neg_control.empty() )
        {
            // Filter to only include negative scores specified by neg_filter
            neg_scores.scores = neg_scores.scores.filter_cols( neg_filter );
        }

    sample_size = original_scores.sample_names.size();

    if( n_opts->approach == "size_factors" )
        {
            std::vector<double> norm_factors( sample_size, 0 );
            compute_size_factors( &norm_factors, &original_scores.scores );
            // normalize the counts
            normalize_counts( &original_scores.scores, &norm_factors );
        }
    else if( n_opts->approach == "neg_diff" )
        {
            matrix<double> neg_norm_diffs;
            compute_neg_diff( &neg_norm_diffs, &neg_scores, &original_scores.scores );
        }
    else if( n_opts->approach == "diff_ratio" )
        {
            matrix<double> neg_norm_diff_ratios;
            compute_neg_diff_ratio( &neg_norm_diff_ratios, &neg_scores, &original_scores.scores );
        }
    else // Column sum default
        {
            std::vector<double> norm_factors( sample_size, 0 );
            get_sum( &norm_factors, &original_scores.scores );
            constant_factor_normalization( &norm_factors, ONE_MILLION );
            // normalize the counts
            normalize_counts( &original_scores.scores, &norm_factors );
        }

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


void module_normalize::constant_factor_normalization( std::vector<double> *cols,
                                                      const std::size_t factor
                                                    )
{
    for( std::size_t index = 0; index < cols->size(); ++index )
        {
            (*cols)[ index ] /= factor;
        }
}

void module_normalize::get_sum( std::vector<double> *dest,
                                matrix<double> *src
                              )
{
    std::size_t index = 0;

    for( index = 0; index < src->nrows(); ++index )
        {
            std::size_t inner_index = 0;
            for( inner_index = 0; inner_index < src->ncols(); ++inner_index )
                {
                    (*dest)[ inner_index ] += (*src)( index, inner_index );
                }

        }
}

void module_normalize::normalize_counts( matrix<double> *original_scores,
                                         const std::vector<double> *norm_factors
                                       )
{
    for( std::size_t index = 0; index < original_scores->nrows(); ++index )
        {
            for( std::size_t inner_index = 0; inner_index < original_scores->ncols(); ++inner_index )
                {
                    (*original_scores)( index, inner_index ) /= (*norm_factors)[ inner_index ];
                }
        }
}

void module_normalize::compute_neg_diff( matrix<double> *norm_diffs,
                                         peptide_score_data_sample_major *neg_scores,
                                         const matrix<double> *data )
{
    double neg_mean = 0;
    // for each sample in data
    for( std::size_t curr_samp = 0; curr_samp < data->ncols(); ++curr_samp )
        {
            // if sb data is given, get mean from current samples rows in neg_scores
            if( !neg_scores->file_name.empty() )
                {
                    neg_mean = stats::geom_mean( data->row_begin( curr_samp ),
                                                data->row_end( curr_samp )
                                                );
                }
            // if not using neg_scores, get mean from current sample rows in input data
                // *will need to filter the input data to only count sample names from neg_filter
            // Now, for each peptide (row) in the sample (column) calculate the difference
            for( std::size_t curr_pep = 0; curr_pep < data->nrows(); ++curr_pep )
                {
                    // store the result in that coordinate of the norm_diffs
                    (*norm_diffs)( curr_pep, curr_samp ) = (*data)(curr_pep,curr_samp) - neg_mean;
                }
        }
}

void module_normalize::compute_neg_diff_ratio( matrix<double> *norm_diff_ratios,
                                         peptide_score_data_sample_major *neg_scores,
                                         const matrix<double> *data )
{
    double neg_mean = 0;
    // for each sample in data
    for( std::size_t curr_samp = 0; curr_samp < data->ncols(); ++curr_samp )
        {
            // if sb data is given, get mean from current samples rows in neg_scores
            if( !neg_scores->file_name.empty() )
                {
                    neg_mean = stats::geom_mean( data->row_begin( curr_samp ),
                                                data->row_end( curr_samp )
                                                );
                }
            // if not using neg_scores, get mean from current sample rows in input data
                // *will need to filter the input data to only count sample names from neg_filter
            // Now, for each peptide (row) in the sample (column) calculate the difference ratio
            for( std::size_t curr_pep = 0; curr_pep < data->nrows(); ++curr_pep )
                {
                    // store the result in that coordinate of the norm_diffs
                    (*norm_diff_ratios)( curr_pep, curr_samp ) =
                                 ((*data)(curr_pep,curr_samp) - neg_mean)/neg_mean;
                }
        }

}

void module_normalize::compute_size_factors( std::vector<double> *size_factors,
                                             const matrix<double> *data
                                           )
{
    std::vector<double> row;
    std::vector<std::vector<double>> geom_mean_factors;
    geom_mean_factors.reserve( data->ncols() );

    std::size_t index = 0;

    for( index = 0; index < data->nrows(); ++index )
        {
            geom_mean_factors.emplace_back( std::vector<double>( data->nrows(), 0 ) ) ;
        }

    for( index = 0; index < data->nrows(); ++index )
        {
            double g_mean = stats::geom_mean( data->row_begin( index ),
                                              data->row_end( index )
                                            );

            for( std::size_t inner_index = 0; inner_index < data->ncols(); ++inner_index )
                {
                    // transpose the data here so computing the medians is easier
                    geom_mean_factors[ inner_index ][ index ] = (*data)( index, inner_index ) / g_mean;
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
            size_factors->push_back( geom_mean_factors[ index ][ median_loc ] );
        }

}

