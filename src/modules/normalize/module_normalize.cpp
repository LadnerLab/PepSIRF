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

#include "logger.h"
#include "module_normalize.h"
#include "options_normalize.h"
#include "time_keep.h"
#include "stats.h"
#define ONE_MILLION 1000000

module_normalize::module_normalize()
{
    name = "Norm";
}

void module_normalize::run( options *opts )
{
    options_normalize *n_opts = (options_normalize*) opts;
    time_keep::timer timer;

    timer.start();

    Log::info("Norm module has started!\n");

    std::string scores_fname = n_opts->peptide_scores_fname;

    omp_set_num_threads( n_opts->num_threads );

    peptide_score_data_sample_major original_scores;
    peptide_score_data_sample_major neg_scores;
    peptide_scoring::parse_peptide_scores( original_scores, scores_fname );
    original_scores.scores = original_scores.scores.transpose();

    std::unordered_set<std::string> neg_filter;
    std::size_t sample_size;

    // Use negative control if provided
    if( !n_opts->neg_control.empty() )
        {
            peptide_scoring::parse_peptide_scores( neg_scores, n_opts->neg_control );
            neg_scores.scores = neg_scores.scores.transpose();
        }
    if(    n_opts->approach == "diff"
        || n_opts->approach == "ratio"
        || n_opts->approach == "diff_ratio" )
        {
            if( !n_opts->neg_names.empty() && !n_opts->neg_id.empty() )
                {
                    std::cout << "Norm does not currently support usage of both negative id [--negative_id,-s] "
                                    "and negative names [--negative_names,-n]. Negative names will be used.\n";
                    boost::split( neg_filter, n_opts->neg_names, boost::is_any_of( "," ) );
                }
            else if( !n_opts->neg_names.empty() )
                {
                    boost::split( neg_filter, n_opts->neg_names, boost::is_any_of( "," ) );
                }
            else if( !n_opts->neg_id.empty() )
                {
                    if( !n_opts->neg_control.empty() )
                        {
                            filter_neg_control_start( &neg_scores, &neg_filter, n_opts->neg_id );
                        }
                    else // if data matrix of sb/neg samples is not provided, use input matrix
                        {
                            filter_neg_control_start( &original_scores, &neg_filter, n_opts->neg_id );
                        }
                }
            else
                {
                    throw std::runtime_error( "ERROR: Must use approach for identifying negative controls. "
                                            "Either a negative id [--negative_id,-s] or negative names [--negative_names,-n].\n" );
                }
        }

    // vector of averages for each peptide to be used in diff, ratio, and diff-ratio approaches.
    std::unordered_map<std::string,double> peptide_averages;
    sample_size = original_scores.sample_names.size();

    if( n_opts->approach == "diff" )
        {
            if( n_opts->neg_control.empty() )
                get_neg_average( &original_scores, &neg_filter, &peptide_averages );
            else
                get_neg_average( &neg_scores, &neg_filter, &peptide_averages );

            compute_diff( &original_scores, &peptide_averages );
        }
    else if( n_opts->approach == "ratio" )
        {
            if( n_opts->neg_control.empty() )
                get_neg_average( &original_scores, &neg_filter, &peptide_averages );
            else
                get_neg_average( &neg_scores, &neg_filter, &peptide_averages );

            compute_ratio( &original_scores, &peptide_averages );
        }
    else if( n_opts->approach == "diff_ratio" )
        {
            if( n_opts->neg_control.empty() )
                get_neg_average( &original_scores, &neg_filter, &peptide_averages );
            else
                get_neg_average( &neg_scores, &neg_filter, &peptide_averages );

            compute_diff_ratio( &original_scores, &peptide_averages );
        }
    else if( n_opts->approach == "size_factors" )
        {
            std::vector<double> norm_factors( sample_size, 0 );
            compute_size_factors( &norm_factors, &original_scores.scores );
            // normalize the counts
            normalize_counts( &original_scores.scores, &norm_factors );
        }
    else if( n_opts->approach == "col_sum" )
        {
            std::vector<double> norm_factors( sample_size, 0 );
            get_sum( &norm_factors, &original_scores.scores );
            constant_factor_normalization( &norm_factors, ONE_MILLION );
            // normalize the counts
            normalize_counts( &original_scores.scores, &norm_factors );
        }
    else
        {
            throw std::runtime_error( "ERROR: Provided approach '" + n_opts->approach + "'. Valid approach must be specified.\n"
            "Available approaches:\n(Default) col_sum\nsize_factors\ndiff\nratio\ndiff_ratio\n"
            "See norm [--help,-h] for further information.\n" );
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

void module_normalize::get_neg_average( peptide_score_data_sample_major *control,
                                        std::unordered_set<std::string> *neg_filter,
                                        std::unordered_map<std::string,double> *pep_averages )
{
    for( std::size_t curr_pep = 0; curr_pep < control->pep_names.size(); ++curr_pep )
        {
            double n_sum = 0.0;
            std::size_t total_samples = 0;
            for( std::size_t curr_samp = 0; curr_samp < control->sample_names.size(); ++curr_samp )
                {
                    if( neg_filter->find( control->sample_names.at( curr_samp ) ) != neg_filter->end() )
                        {
                            n_sum += control->scores.at( curr_pep, curr_samp );
                            ++total_samples;
                        }
                }
            pep_averages->emplace(control->pep_names.at(curr_pep), n_sum/total_samples);
        }
}

void module_normalize::filter_neg_control_start( peptide_score_data_sample_major *score_data,
                                                 std::unordered_set<std::string> *neg_filter,
                                                 std::string start )
{
    for( const auto& sample : score_data->sample_names )
        {
            if( sample.rfind( start, start.length() ) != std::string::npos )
                {
                    neg_filter->emplace( sample );
                }
        }
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

void module_normalize::compute_diff( peptide_score_data_sample_major *norm_diffs,
                                     std::unordered_map<std::string,double> *peptide_avgs )
{
    for( const auto& name_score_pair : *peptide_avgs )
        {
            std::string peptide_name = name_score_pair.first;
            double neg_mean = name_score_pair.second;
            for( std::size_t curr_samp = 0; curr_samp < norm_diffs->scores.ncols(); ++curr_samp )
                {
                    // store the result in that coordinate of the norm_diffs
                    norm_diffs->scores.at( peptide_name, norm_diffs->sample_names.at(curr_samp) ) = 
                    norm_diffs->scores.at( peptide_name, norm_diffs->sample_names.at(curr_samp) ) - neg_mean;
                                        
                }
        }
}

void module_normalize::compute_ratio( peptide_score_data_sample_major *norm_ratios,
                                     std::unordered_map<std::string,double> *peptide_avgs )
{
    for( const auto& name_score_pair : *peptide_avgs )
        {
            std::string peptide_name = name_score_pair.first;
            double neg_mean = name_score_pair.second;
            for( std::size_t curr_samp = 0; curr_samp < norm_ratios->scores.ncols(); ++curr_samp )
                {
                    // store the result in that coordinate of the norm_diffs
                    norm_ratios->scores.at( peptide_name, norm_ratios->sample_names.at(curr_samp) ) =
                    norm_ratios->scores.at( peptide_name, norm_ratios->sample_names.at(curr_samp) )/neg_mean;
                }
        }
}

void module_normalize::compute_diff_ratio( peptide_score_data_sample_major *norm_diff_ratios,
                                           std::unordered_map<std::string,double> *peptide_avgs )
{
    for( const auto& name_score_pair : *peptide_avgs )
        {
            std::string peptide_name = name_score_pair.first;
            double neg_mean = name_score_pair.second;
            for( std::size_t curr_samp = 0; curr_samp < norm_diff_ratios->scores.ncols(); ++curr_samp )
                {
                    // store the result in that coordinate of the norm_diffs
                    norm_diff_ratios->scores.at( peptide_name, norm_diff_ratios->sample_names.at(curr_samp) ) = 
                    (norm_diff_ratios->scores.at( peptide_name, norm_diff_ratios->sample_names.at(curr_samp) ) - neg_mean)/neg_mean;
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

