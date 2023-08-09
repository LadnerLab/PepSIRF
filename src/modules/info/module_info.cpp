#include <iostream>
#include <fstream>
#include <numeric>
#include <iomanip>
#include <boost/algorithm/find_backward.hpp> //TODO: include in header file

#include "logger.h"
#include "module_info.h"
#include "peptide_scoring.h"
#include "file_io.h"

module_info::module_info()
{
    name = "Info";
}

void module_info::run( options *opts )
{

    options_info &i_opts = *(options_info*) opts;
    peptide_score_data_sample_major scores;

    peptide_scoring::parse_peptide_scores( scores,
                                           i_opts.in_fname
                                         );

    Log::info(
        "Number of samples: "
        + std::to_string(scores.sample_names.size()) + "\n"
    );
    Log::info(
        "Number of peptides: "
        + std::to_string(scores.pep_names.size()) + "\n"
    );

    if( !i_opts.out_samples_fname.empty() )
        {
            std::ofstream sample_names{ i_opts.out_samples_fname };

            pepsirf_io::write_file( sample_names,
                                    scores.sample_names.begin(),
                                    scores.sample_names.end(),
                                    "\n"
                                   );
        }
    if( !i_opts.out_pep_names_fname.empty() )
        {
            std::ofstream peptide_names{ i_opts.out_pep_names_fname };

            pepsirf_io::write_file( peptide_names,
                                    scores.pep_names.begin(),
                                    scores.pep_names.end(),
                                    "\n"
                                   );
        }

    if( !i_opts.out_col_sums_fname.empty() )
        {
            std::ofstream peptide_names{ i_opts.out_col_sums_fname };
            peptide_names << "Sample name\tSum of probe scores\n";

            auto& matr = scores.scores;

            for( const auto &sample_n : matr.get_row_labels() )
                {
                    double sum = std::accumulate( matr.row_begin( sample_n.second ),
                                                  matr.row_end( sample_n.second ),
                                                  static_cast<double>( 0.0 )
                                                );

                    peptide_names << sample_n.first
                                  << "\t"
                                  << std::fixed << std::setprecision( 2 ) << sum
                                  << "\n";
                }

        }

    if ( !i_opts.in_replicates_fname.empty() && !i_opts.out_avgs_fname.empty() )
        {
            std::ifstream replicate_names{ i_opts.in_replicates_fname };
            std::ofstream averages{ i_opts.out_avgs_fname };

            std::unordered_map<std::string, std::vector<std::string>> name_file_samples = {};

            // Create map of sample names with the base sample name as the key
            // And a list of its associated sample names as the value
            for ( std::string current_line; std::getline(replicate_names, current_line); )
                {
                    // Get the base sample name and remove it from the current line
                    std::string delimiter = "\t";
                    size_t pos = current_line.find( delimiter );
                    std::string base_sample = current_line.substr( 0, pos );
                    current_line.erase( 0, pos + delimiter.length() );

                    // Loop through the rest of the sample names on the line
                    // And add them to a vector
                    std::vector<std::string> sample_list;
                    while( ( pos = current_line.find( delimiter ) ) != std::string::npos )
                        {
                            std::string sample = current_line.substr( 0, pos );
                            sample_list.emplace_back( sample );
                            current_line.erase( 0, pos + delimiter.length() );
                        }
                    // Add the last sample to the vector
                    sample_list.emplace_back( current_line.substr( 0, current_line.find( "\r" ) ) );
                    name_file_samples.emplace( std::make_pair( base_sample, sample_list ) );
                }

            // Check that sample names in names file match those in the input file
            bool invalid_sample_found;
            std::vector<std::string> invalid_samples = {};
            std::vector<std::string> found_samples = {};
            for ( std::size_t i = 0; i < scores.sample_names.size(); i++ )
                {
                    invalid_sample_found = true;

                    // Loop through samples found in name file
                    for ( auto samples : name_file_samples )
                        {
                            // loop through replicate samples provided in sample
                            for ( std::string sample : samples.second )
                                {
                                    // If samples in input & name files match, write to output file
                                    if ( std::find(
                                            scores.sample_names.begin(),
                                            scores.sample_names.end(),
                                            sample
                                            ) != scores.sample_names.end()
                                        )
                                        {
                                            invalid_sample_found = false;
                                            found_samples.emplace_back( sample );
                                        }
                                    // otherwise, check sample was not added to invalid samples
                                    else if ( boost::algorithm::find_backward(
                                                invalid_samples.begin(),
                                                invalid_samples.end(),
                                                samples.first
                                                ) == invalid_samples.end()
                                            )
                                        {
                                            // add sample name to invalid samples and write warning
                                            invalid_samples.emplace_back( samples.first );
                                            Log::warn("Invalid replicate name found in " + samples.first + "\n");
                                        }
                                }
                        }
                }

            // verify no duplicate replicates
            bool duplicate_samples_found = false;
            std::vector<std::string> duplicate_samples;
            found_samples = {};
            // loop through replicates in input matrix
            for ( auto sample : scores.sample_names )
                {
                    // check replicate already found
                    if ( std::find(
                            found_samples.begin(),
                            found_samples.end(),
                            sample
                            ) != found_samples.end()
                        && boost::algorithm::find_backward(
                            duplicate_samples.begin(),
                            duplicate_samples.end(),
                            sample
                            ) == duplicate_samples.end()
                        )
                        {
                            duplicate_samples_found = true;
                            duplicate_samples.emplace_back( sample );
                        }
                    // otherwise, assume unique replicate name
                    else
                        {
                            found_samples.emplace_back( sample );
                        }
                }

            // Print warning with all duplicate samples found
            if ( duplicate_samples_found )
                {
                    Log::warn("Duplicate samples in input file:\n");

                    for ( std::string sample : duplicate_samples )
                        {
                            Log::info(sample + "\n");
                        }
                }

            // Write base sample names as row headers in output file
            averages << "Sequence name";
            for( auto sample : name_file_samples )
                {
                    if ( std::find( invalid_samples.begin(), invalid_samples.end(), sample.first ) == invalid_samples.end() )
                        {
                            averages << "\t" << sample.first;
                        }
                }
            averages << "\n";

            // calculate sample averages and write to output matrix
            double rep_total = 0.0;
            // loop over peptides from input matrix
            for ( size_t pep_idx = 0;
                  pep_idx < scores.pep_names.size();
                  pep_idx += 1
                )
                {
                    averages << scores.pep_names[ pep_idx ];

                    // loop over provided replicates
                    for ( auto sample : name_file_samples )
                        {
                            rep_total = 0.0;

                            // loop over replicate from input matrix
                            for ( size_t rep_idx = 0;
                                  rep_idx < scores.sample_names.size();
                                  rep_idx += 1
                                )
                                {
                                    // check current replicate included in sample
                                    if ( std::find(
                                            sample.second.begin(),
                                            sample.second.end(),
                                            scores.sample_names[ rep_idx ]
                                            ) == sample.second.end()
                                       )
                                        {
                                            continue;  // move to next replicate from input matrix
                                        }

                                    // add scores to replicate total
                                    rep_total += scores.scores.at( rep_idx, pep_idx );
                                }

                            // average replicate total and print to output matrix
                            averages << "\t" << std::fixed << std::setprecision( 6 )
                                     << rep_total / (double)( sample.second.size() );
                        }

                    averages << "\n";
                }
        }
}

