#include <boost/algorithm/string.hpp>
#include <stdexcept>
#include <numeric>
#include <fstream>
#include <iostream>

#include "module_p_enrich.h"
#include "time_keep.h"
#include "omp_opt.h"
#include "fs_tools.h"
#include "predicate.h"
#include "file_io.h"
#include "setops.h"


void module_p_enrich::run( options *opts )
{
    options_p_enrich *p_opts = (options_p_enrich*) opts;
    time_keep::timer timer;
    timer.start();

    peptide_score_data_sample_major raw_scores;
    peptide_score_data_sample_major *raw_scores_ptr = nullptr;
    std::vector<double> raw_score_params;
    std::vector<std::pair<peptide_score_data_sample_major,std::vector<double>>> matrix_thresh_pairs;
    // storing all scores in from each matrix file provided in matrix_scores
    matrix_thresh_pairs.resize( p_opts->matrix_thresh_fname_pairs.size() );

    std::size_t curr_matrix;

    for( curr_matrix = 0; curr_matrix < p_opts->matrix_thresh_fname_pairs.size(); curr_matrix++ )
        {
            matrix_thresh_pairs[curr_matrix] = std::make_pair( peptide_score_data_sample_major(), std::vector<double>{ 0.0, 0.0 } );
            peptide_scoring::parse_peptide_scores( matrix_thresh_pairs[curr_matrix].first, p_opts->matrix_thresh_fname_pairs[ curr_matrix ].first );
            std::vector<std::string> str_thresholds;
            boost::split( str_thresholds, p_opts->matrix_thresh_fname_pairs[ curr_matrix ].second, boost::is_any_of( "," ) );
            std::transform( str_thresholds.begin(), str_thresholds.end(), matrix_thresh_pairs[curr_matrix].second.begin(),
                [&]( const std::string& val)
                {
                    return std::stod( val );
                });
        }

    std::ifstream pairs_file;
    pairs_file.exceptions( std::ios::badbit );

    try
        {
            pairs_file.open( p_opts->in_samples_fname );
        }

    catch( const std::ifstream::failure& e )
        {
            throw std::runtime_error( "Unable to open the provided samples file" );
        }

    auto sample_pairs = parse_samples( pairs_file );

    pairs_file.close();

    bool raw_counts_included = !p_opts->in_raw_scores_fname.empty();

    if( raw_counts_included )
        {
            peptide_scoring::parse_peptide_scores( raw_scores, p_opts->in_raw_scores_fname );
            raw_scores_ptr = &raw_scores;

            std::vector<std::string> str_constraints;
            boost::split( str_constraints, p_opts->raw_scores_params_str, boost::is_any_of( "," ) );
            raw_score_params.resize( str_constraints.size() );
            std::transform( str_constraints.begin(), str_constraints.end(), raw_score_params.begin(),
                [&]( const std::string& val)
                {
                    return std::stod( val );
                });

        }

    auto output_path = fs_tools::path( p_opts->out_dirname );
    bool dir_exists = !fs_tools::create_directories( output_path );

    if( dir_exists )
        {
            std::cout << "WARNING: the directory '" << p_opts->out_dirname
                      << "' exists, any files with "
                      << "colliding filenames will be overwritten!\n";
        }

    using sample_name_set = std::unordered_set<std::string>;
    sample_name_set raw_score_sample_names{ raw_scores.sample_names.begin(),
            raw_scores.sample_names.end()
            };


    bool valid_candidate;
    bool raw_count_enriched = true;
    peptide_score_data_sample_major *curr_matrix_ptr;
    for( std::size_t sample_idx = 0; sample_idx < sample_pairs.size(); ++sample_idx )
        {
            std::vector<paired_score> enriched_probes;
            std::vector<std::map<std::string,paired_score>> all_enrichment_candidates;
            all_enrichment_candidates.resize( matrix_thresh_pairs.size() );
            for( curr_matrix = 0; curr_matrix < matrix_thresh_pairs.size(); ++curr_matrix )
                {

                    curr_matrix_ptr = &matrix_thresh_pairs[curr_matrix].first;

                    // list of samples
                    sample_name_set matrix_sample_names{ curr_matrix_ptr->sample_names.begin(),
                        curr_matrix_ptr->sample_names.end() };

                    // If included, raw_count_sample_names must also equal
                    //  names, otherwise we do not care.

                    if( raw_counts_included &&
                    ( raw_score_sample_names.size() != matrix_sample_names.size() ) )
                        {
                            throw std::runtime_error( "The samplenames provided in each input file are "
                                                    "not the same"
                                                    );
                        }
                    bool first_in = matrix_sample_names.find( sample_pairs[sample_idx].first )
                                != matrix_sample_names.end();
                    bool second_in = matrix_sample_names.find( sample_pairs[sample_idx].second )
                                != matrix_sample_names.end();

                    if( first_in ^ second_in )
                        {
                            std::string not_found = first_in ? sample_pairs[sample_idx].second
                                : sample_pairs[sample_idx].first;
                            std::string error_msg = "The sample '"
                                + not_found + "'"
                                + " was not found in any of the matrix files. ";

                            throw std::runtime_error( error_msg );
                        }
                    else if( !( first_in && second_in ) )
                        {
                            std::cout << "WARNING: the samples '"
                                << sample_pairs[sample_idx].first
                                << "' and '"
                                << sample_pairs[sample_idx].second
                                << "' were not in '"
                                << curr_matrix_ptr->file_name
                                << "' and will be skipped. \n";
                        }
                    else
                        {
                            //create a vector of enrichment candidates - each vector of paired scores for a thresholds file provided
                            get_enrichment_candidates( &all_enrichment_candidates[curr_matrix],
                                                       curr_matrix_ptr,
                                                       raw_scores_ptr,
                                                       sample_pairs[ sample_idx ]
                                                     );
                            // only attempt to get the column sum if raw counts
                            // have been specified
                            std::pair<double,double>
                                col_sums{ 0.0, 0.0 };

                            if( raw_counts_included )
                                {
                                    col_sums = get_raw_sums( all_enrichment_candidates[curr_matrix].begin(),
                                                    all_enrichment_candidates[curr_matrix].end()
                                                );
                                }
                            // determine if raw count is enriched if meets or is over threshold
                            raw_count_enriched = raw_counts_included
                                ? ( raw_count_enriched && thresholds_met( col_sums, raw_score_params ) )
                                : true;
                        }
                }
            if( raw_count_enriched )
                {
                    // verify all candidates of the same peptide name are over given thresholds
                    for( const auto& candidate : all_enrichment_candidates[0] )
                        {
                            std::string pep_name = candidate.first;
                            valid_candidate = true;
                            for( std::size_t curr_map = 0; curr_map < all_enrichment_candidates.size(); ++curr_map )
                                {
                                    // if( all_enrichment_candidates[curr_map].find( pep_name ) != all_enrichment_candidates[curr_map].end() )
                                    //     {

                                    if( valid_candidate
                                        && thresholds_met( all_enrichment_candidates[curr_map].at( pep_name ).score,
                                                        matrix_thresh_pairs[curr_map].second ) )
                                        {
                                            valid_candidate = true;
                                        }
                                    else
                                        {
                                            valid_candidate = false;
                                        }
                                        // }
                                }
                            if( valid_candidate )
                                {
                                    paired_score enriched_probe = { pep_name,
                                                                    candidate.second.score,
                                                                    candidate.second.raw_score };
                                    enriched_probes.emplace_back( enriched_probe );
                                }
                        }

                }
            if( !enriched_probes.empty() )
                {
                    std::string outf_name =
                        p_opts->out_dirname + '/'
                        + sample_pairs[ sample_idx ].first
                        + p_opts->out_fname_join
                        + sample_pairs[ sample_idx ].second
                        + p_opts->out_suffix;
                    std::ofstream out_file{ outf_name, std::ios_base::out };

                    pepsirf_io::write_file( out_file,
                                            enriched_probes.begin(),
                                            enriched_probes.end(),
                                            "\n"
                                        );
                }

        }
}

std::vector<module_p_enrich::sample_type>
module_p_enrich::parse_samples( std::istream& file )
{
    std::vector<sample_type> return_val;

    std::string current_line;
    std::vector<std::string> values;
    // reserve the expected amount
    values.reserve( 2 );

    while( std::getline( file, current_line ) )
        {
            if( current_line.find("\r") != std::string::npos )
                current_line.erase( current_line.find( "\r" ) );
            boost::split( values,
                          current_line,
                          boost::is_any_of( "\t" )
                       );

            if( values.size() != 2 )
                {
                    throw std::runtime_error( "The input samples file must "
                                              "contain two samples per line "
                                            );
                }

            return_val.emplace_back( std::make_pair( values[ 0 ],
                                                     values[ 1 ]
                                                   )
                                   );
        }

    return return_val;
}

void module_p_enrich::get_enrichment_candidates( std::map<std::string,paired_score> *enrichment_candidates,
                                            const peptide_score_data_sample_major *matrix_score_data,
                                            const peptide_score_data_sample_major *raw_score_data,
                                            const std::pair<std::string,std::string> sample_names
                                          )
{
    std::vector<paired_score> candidates;
    std::size_t pep_idx = 0;

    using pair = std::pair<double,double>;

    for( pep_idx = 0; pep_idx < matrix_score_data->pep_names.size(); ++pep_idx )
        {
            std::string pep_name = matrix_score_data->pep_names[ pep_idx ];
            // get matrix score (norm or zscore) in each of the sample
            pair matrix_scores = { matrix_score_data->scores( sample_names.first, pep_name ),
                             matrix_score_data->scores( sample_names.second, pep_name )
                           };

            pair raw_scores = { 0.0, 0.0 };

            if( raw_score_data != nullptr)
                {
                    raw_scores = { raw_score_data
                                   ->scores( sample_names.first, pep_name ),
                                   raw_score_data
                                   ->scores( sample_names.second, pep_name )
                                 };
                }

            paired_score candidate{ matrix_scores,
                                    raw_scores
                                  };

            enrichment_candidates->emplace( pep_name, candidate );
        }

}
