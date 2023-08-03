#include <boost/algorithm/string.hpp>
#include <stdexcept>
#include <numeric>
#include <fstream>
#include <iostream>
#include <string>

#include "logger.h"
#include "module_enrich.h"
#include "time_keep.h"
#include "omp_opt.h"
#include "fs_tools.h"
#include "predicate.h"
#include "file_io.h"
#include "setops.h"

module_enrich::module_enrich()
{
    name = "Enrich";
}

void module_enrich::run( options *opts )
{
    options_enrich *e_opts = (options_enrich*) opts;
    time_keep::timer timer;
    timer.start();

    peptide_score_data_sample_major raw_scores;
    peptide_score_data_sample_major *raw_scores_ptr = nullptr;
    std::vector<double> raw_score_params;
    std::vector<std::pair<peptide_score_data_sample_major,std::vector<double>>> matrix_thresh_pairs;
    std::unordered_map<std::string,std::string> enrichment_failures;
    matrix_thresh_pairs.resize( e_opts->matrix_thresh_fname_pairs.size() );
    bool shorthand_output_filenames = e_opts->truncate_names;
    bool low_reads = e_opts->low_raw_reads;
    std::size_t curr_matrix;

    for( curr_matrix = 0; curr_matrix < e_opts->matrix_thresh_fname_pairs.size(); curr_matrix++ )
        {
            matrix_thresh_pairs[curr_matrix] = std::make_pair( peptide_score_data_sample_major(), std::vector<double>{ 0.0 } );
            peptide_scoring::parse_peptide_scores( matrix_thresh_pairs[curr_matrix].first, e_opts->matrix_thresh_fname_pairs[ curr_matrix ].first );
            std::vector<std::string> str_thresholds;
            boost::split( str_thresholds, e_opts->matrix_thresh_fname_pairs[ curr_matrix ].second, boost::is_any_of( "," ) );
            matrix_thresh_pairs[curr_matrix].second.resize( str_thresholds.size() );
            std::transform( str_thresholds.begin(), str_thresholds.end(), matrix_thresh_pairs[curr_matrix].second.begin(),
                [&]( const std::string& val)
                {
                    return std::stod( val );
                });
        }

    std::ifstream samples_file;
    std::vector<module_enrich::sample_type> samples_list;
    if( !e_opts->in_samples_fname.empty() )
        {
            samples_file.open( e_opts->in_samples_fname );
            if( samples_file.fail() )
                {
                    Log::error("Unable to open the provided samples file!");
                }

            samples_list = parse_samples( samples_file );

            samples_file.close();

        }
    // if no samples file is given, then all samples will be included and treated as samples assayed in single replicate.
    if( samples_list.empty() )
        {
            std::for_each( matrix_thresh_pairs[0].first.sample_names.begin(), matrix_thresh_pairs[0].first.sample_names.end(),
                    [&]( std::string sample_name )
                        {
                            samples_list.emplace_back( module_enrich::sample_type{sample_name} );
                        });
        }
    auto output_path = fs_tools::path( e_opts->out_dirname );
    bool dir_exists = !fs_tools::create_directories( output_path );

    if( dir_exists )
        {
            Log::warn(
                "The directory '" + e_opts->out_dirname + "' exists, any files"
                + " with colliding filenames will be overwritten!\n"
            );
        }

    bool raw_counts_included = !e_opts->in_raw_scores_fname.empty();

    if( raw_counts_included )
        {
            peptide_scoring::parse_peptide_scores( raw_scores, e_opts->in_raw_scores_fname );
            raw_scores_ptr = &raw_scores;

            std::vector<std::string> str_constraints;
            boost::split( str_constraints, e_opts->raw_scores_params_str, boost::is_any_of( "," ) );
            raw_score_params = {0.0};
            std::transform( str_constraints.begin(), str_constraints.end(), raw_score_params.begin(),
                [&]( const std::string& val)
                {
                    return std::stod( val );
                });
        }
    using sample_name_set = std::set<std::string>;
    sample_name_set raw_score_sample_names{ raw_scores.sample_names.begin(),
                                            raw_scores.sample_names.end()
                                            };

    std::vector<std::string> problem_replicates;
    std::vector<std::vector<std::string>> removed_reps;
    peptide_score_data_sample_major *curr_matrix_ptr;
    for( std::size_t sample_idx = 0; sample_idx < samples_list.size(); ++sample_idx )
        {
            bool raw_count_enriched = true;
            std::vector<std::string> enriched_probes;

            // raw_score_lists: a list of peptide scores, and their sample
            // indexed with: x-axis is the peptide, y-axis is the sample name
            std::vector<std::vector<double>> raw_score_lists;

            std::vector<std::map<std::string,std::vector<double>>> all_enrichment_candidates;
            all_enrichment_candidates.resize( matrix_thresh_pairs.size() );
            raw_score_lists.resize( raw_scores.scores.ncols() );
            std::vector<double> col_sums;

            if( raw_counts_included )
                {
                    get_raw_scores( &raw_score_lists, raw_scores_ptr, samples_list[sample_idx] );
                    col_sums = get_raw_sums( raw_score_lists );
                }

            if( low_reads )
                {
                    std::vector<std::string> saved_reps;

                    std::vector<double>::iterator col_sum = col_sums.begin();
                    for ( std::vector<std::string>::iterator rep = samples_list[ sample_idx ].begin();
                            rep != samples_list[ sample_idx ].end(); )
                        {
                            if ( *col_sum < *std::min_element( raw_score_params.begin(), raw_score_params.end() ) )
                                {
                                    saved_reps.emplace_back( *rep );
                                    rep = samples_list[ sample_idx ].erase( rep );
                                    col_sum = col_sums.erase( col_sum );
                                }
                            else
                                {
                                    rep += 1;
                                    col_sum += 1;
                                }
                        }

                    if ( !saved_reps.empty() ) // keeps from having to account when outputing to enrich failure file
                        {
                            removed_reps.emplace_back( saved_reps );
                        }
                }

            raw_count_enriched = raw_counts_included
                ? ( !col_sums.empty() && raw_count_enriched && thresholds_met( col_sums, raw_score_params ) )
                : true;

            if( raw_count_enriched )
            {
                // Get candidates for each input matrix for the current sample
                for( curr_matrix = 0; curr_matrix < matrix_thresh_pairs.size(); ++curr_matrix )
                    {
                        curr_matrix_ptr = &matrix_thresh_pairs[curr_matrix].first;

                        sample_name_set matrix_sample_names{ curr_matrix_ptr->sample_names.begin(),
                            curr_matrix_ptr->sample_names.end() };
                        sample_name_set samples_sample_names{ samples_list[sample_idx].begin(), samples_list[sample_idx].end() };

                        if( raw_counts_included &&
                        ( raw_score_sample_names.size() != matrix_sample_names.size() ) )
                            {
                                Log::error(
                                    // y'all still wanna defend C++ as a good lang?
                                    std::string("The samplenames provided in each input")
                                    + " file are not the same!\n"
                                );
                            }
                        std::vector<std::string> sample_diffs;
                        std::set_difference( samples_sample_names.begin(), samples_sample_names.end(),
                                             matrix_sample_names.begin(), matrix_sample_names.end(), sample_diffs.begin() );
                        // if some samples are not found that are provided by samples option, then throw an error
                        if( !sample_diffs.empty() && sample_diffs.size() != samples_sample_names.size() )
                            {
                                Log::info(
                                    std::string("The listed samples were not found in")
                                    + " matrix '" + curr_matrix_ptr->file_name
                                    + "'.\n"
                                );

                                for( auto sample = sample_diffs.begin(); sample != sample_diffs.end(); sample++ )
                                    {
                                        // TODO: send through Log::info()
                                        std::cout << *sample << "\n";
                                    }

                                Log::error(
                                    std::string("Verify the correct sample names are")
                                    + " provided in (--samples, -s).\n"
                                );

                            }
                        // otherwise, if all of the samples are not found that are
                        // provided by samples option, then give a warning
                        // that these samples were not in the matrix.
                        else if( sample_diffs.size() == samples_sample_names.size() )
                            {
                                Log::warn(
                                    std::string("The listed samples were not found in")
                                    + " matrix '" + curr_matrix_ptr->file_name
                                    + "' and will be skipped.\n"
                                );

                                for( auto sample = sample_diffs.begin(); sample != sample_diffs.end(); sample++ )
                                    {
                                        // TODO: send through Log::info()
                                        std::cout << *sample << "\n";
                                    }
                            }
                        else
                            {
                                get_enrichment_candidates(
                                    &all_enrichment_candidates[curr_matrix],
                                    curr_matrix_ptr,
                                    samples_list[ sample_idx ]
                                );
                            }
                    }

                for( const auto& candidate : all_enrichment_candidates[0] )
                    {
                        std::string pep_name = candidate.first;
                        std::size_t valid_candidates = 0;
                        std::vector<std::map<std::string, std::vector<double>>> ret;
                        for( std::size_t curr_map = 0; curr_map < all_enrichment_candidates.size(); ++curr_map )
                            {
                                if (thresholds_met(
                                        all_enrichment_candidates[curr_map].at(pep_name),
                                        matrix_thresh_pairs[curr_map].second
                                ))
                                ++valid_candidates;
                            }

                        if( valid_candidates == matrix_thresh_pairs.size() )
                            {
                                enriched_probes.emplace_back( pep_name );
                            }
                    }
            }
            
			if (
				!samples_list[sample_idx].empty()
                && raw_counts_included
                && !raw_count_enriched
                && !e_opts->out_enrichment_failure.empty()
            )
                {
                    std::vector<double>::iterator min_thresh = std::min_element( raw_score_params.begin(), raw_score_params.end() );
                    std::vector<double>::iterator max_thresh = std::max_element( raw_score_params.begin(), raw_score_params.end() );
                    std::vector<double>::iterator min_sum = std::min_element( col_sums.begin(), col_sums.end() );
                    std::vector<double>::iterator max_sum = std::max_element( col_sums.begin(), col_sums.end() );
        
                    std::string samplenames = "";
                    std::string prob_reps_report = "";
                    std::vector<std::string> replicate_names;
                    std::for_each( samples_list[sample_idx].begin(), samples_list[sample_idx].end() - 1,
                        [&]( std::string name )
                            {
                                replicate_names.emplace_back( name );
                                samplenames.append( name + ", " );
                            });
                    replicate_names.emplace_back( *( samples_list[ sample_idx ].end() - 1 ) );
                    samplenames.append( *(samples_list[sample_idx].end() - 1) );
                    enrichment_failures.emplace( samplenames, "raw" );

                    std::string separate = "";
                    for ( std::size_t name_idx = 0; name_idx < replicate_names.size() - 1; name_idx += 1 )
                        {
                            if ( col_sums[ name_idx ] == *min_sum && col_sums[ name_idx ] < *min_thresh )
                                {
                                    prob_reps_report.append( replicate_names[ name_idx ] + " (min)" + separate );
                                }
                            else if ( col_sums[ name_idx ] == *max_sum && col_sums[ name_idx ] < *max_thresh )
                                {
                                    prob_reps_report.append( replicate_names[ name_idx ] + " (max)" + separate );
                                }

                            separate = ", ";
                        }

                    if ( col_sums[ replicate_names.size() - 1 ] == *min_sum
                            && col_sums[ replicate_names.size() - 1 ] < *min_thresh )
                        {
                            if ( !prob_reps_report.empty() )
                                {
                                    prob_reps_report.append( separate + replicate_names[ replicate_names.size() - 1 ] + " (min)" );
                                }
                            else
                                {
                                    prob_reps_report.append( replicate_names[ replicate_names.size() - 1 ] + " (min)" );
                                }
                        }
                    else if ( col_sums[ replicate_names.size() - 1 ] == *max_sum
                                && col_sums[ replicate_names.size() - 1 ] < *max_thresh )
                        {
                            if ( !prob_reps_report.empty() )
                                {
                                    prob_reps_report.append( separate + replicate_names[ replicate_names.size() - 1 ] + " (max)" );
                                }
                            else
                                {
                                    prob_reps_report.append( replicate_names[ replicate_names.size() - 1 ] + " (max)" );
                                }
                        }

                    problem_replicates.emplace_back( prob_reps_report );
                }
            else if (
				!samples_list[sample_idx].empty()
                && enriched_probes.empty()
				&& !low_reads
                && !e_opts->out_enrichment_failure.empty()
			)
                {	// collects replicates with no enriched peptides
                    std::string samplenames = "";
                    std::for_each( samples_list[sample_idx].begin(), samples_list[sample_idx].end() - 1,
                        [&]( std::string name )
                            {
                                samplenames.append( name + ", " );
                            });
                    samplenames.append( *(samples_list[sample_idx].end() - 1) );
                    enrichment_failures.emplace( samplenames, "peptides" );
                }

            std::string outf_name = e_opts->out_dirname + '/';
            if ( shorthand_output_filenames && samples_list[sample_idx].size() > 3 )
                {
                    for( std::size_t name_idx = 0; name_idx < 3; name_idx++ )
                        {
                            outf_name += samples_list[sample_idx][name_idx] + e_opts->out_fname_join;
                        }
                    outf_name += std::to_string(samples_list[sample_idx].size() - 3) + "more" + e_opts->out_suffix;
                }
            else if ( !samples_list[sample_idx].empty() )
                {
                    for( std::size_t name_idx = 0; name_idx < samples_list[sample_idx].size() - 1; name_idx++ )
                        {
                            outf_name += samples_list[sample_idx][name_idx] + e_opts->out_fname_join;
                        }
                    outf_name += *(samples_list[sample_idx].end() - 1) + e_opts->out_suffix;
                }
            std::ofstream out_file{ outf_name, std::ios_base::out };

            pepsirf_io::write_file(
                out_file,
                enriched_probes.begin(),
                enriched_probes.end(),
                "\n"
            );

            if ( enriched_probes.size() == 0 )
                {
                   out_file << ' ';
                }
        }

    if ( !e_opts->out_enrichment_failure.empty() )
        {
            if ( !enrichment_failures.empty() )
                {
                    std::size_t pr_idx = problem_replicates.size() - 1;
                    std::string outf_name = e_opts->out_dirname + '/' + e_opts->out_enrichment_failure;

                    // write to file
                    std::ofstream out_file{ outf_name, std::ios_base::out };
                    out_file << "Replicates\tReason\tProblem Replicates\n";

                    for( auto& line : enrichment_failures )
                        {
                            out_file << line.first;
                            if( line.second == "raw" )
                                {
                                    out_file << "\tRaw read count threshold\t";
                                    out_file << problem_replicates[ pr_idx ] << std::endl;
                                    pr_idx -= 1;
                                }
                            else
                                {
                                    out_file << "\tNo enriched peptides\n";
                                }
                        }

                    out_file.close();
                }
            else if ( !removed_reps.empty() )
                {
                    std::string outf_name = e_opts->out_dirname + '/' + e_opts->out_enrichment_failure;

                    // write to file
                    std::ofstream out_file{ outf_name, std::ios_base::out };
                    out_file << "Removed Replicates\n";

                    for ( std::size_t s = 0; s < removed_reps.size(); s += 1 )
                        {
                            for ( std::size_t r = 0; r < removed_reps[ s ].size() - 1; r += 1 )
                                {
                                    out_file << removed_reps[ s ][ r ] << ", ";
                                }

                            out_file << removed_reps[ s ][ removed_reps[ s ].size() - 1 ] << std::endl;
                        }

                    out_file.close();
                }
        }
}

std::vector<module_enrich::sample_type> module_enrich::parse_samples( std::istream& file )
{
    std::vector<sample_type> return_val;

    std::string current_line;
    std::vector<std::string> values;
    // reserve the min expected amount
    values.reserve( 1 );

    while( std::getline( file, current_line ) )
        {
            if( current_line.find("\r") != std::string::npos )
                current_line.erase( current_line.find( "\r" ) );
            boost::split( values,
                          current_line,
                          boost::is_any_of( "\t" )
                       );
            return_val.emplace_back(values);

        }

    return return_val;
}

void module_enrich::get_enrichment_candidates( std::map<std::string,std::vector<double>> *enrichment_candidates,
                                            const peptide_score_data_sample_major *matrix_score_data,
                                            const std::vector<std::string> sample_names
                                          )
{
    std::vector<paired_score> candidates;
    std::size_t pep_idx;

    for( pep_idx = 0; pep_idx < matrix_score_data->pep_names.size(); ++pep_idx )
        {
            std::string pep_name = matrix_score_data->pep_names[ pep_idx ];
            std::vector<double> matrix_scores;
            for( const auto& name : sample_names )
                {
                    matrix_scores.emplace_back( matrix_score_data->scores( name, pep_name ) );
                }
            enrichment_candidates->emplace( pep_name, matrix_scores );
        }

}

void module_enrich::get_raw_scores( std::vector<std::vector<double>> *raw_scores_dest,
                                      const peptide_score_data_sample_major *raw_score_data,
                                      const std::vector<std::string> sample_names )
{
    for( std::size_t pep_idx = 0; pep_idx < raw_score_data->pep_names.size(); ++pep_idx )
    {
        std::vector<double> raw_scores;

        for( const auto& name : sample_names )
            {
                raw_scores.emplace_back( raw_score_data->scores( name, raw_score_data->pep_names[pep_idx] ) );
            }
        raw_scores_dest->at( pep_idx ) = raw_scores;
    }
}

