#include <iostream>
#include <fstream>
#include <numeric>
#include <iomanip>
#include <boost/algorithm/find_backward.hpp> //TODO: include in header file

#include "module_info.h"
#include "peptide_scoring.h"
#include "file_io.h"

void module_info::run( options *opts )
{

    options_info &i_opts = *(options_info*) opts;
    peptide_score_data_sample_major scores;

    peptide_scoring::parse_peptide_scores( scores,
                                           i_opts.in_fname
                                         );


    std::cout << "Number of samples:  " << scores.sample_names.size() << "\n";
    std::cout << "Number of peptides: " << scores.pep_names.size() << "\n";

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

            for( const auto sample_n : matr.get_row_labels() )
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

    if (!i_opts.in_replicates_fname.empty() && !i_opts.out_avgs_fname.empty())
        {
            std::unordered_map<std::string, std::vector<int>> sample_map = {};
            std::unordered_map<std::string, std::vector<std::string>> name_file_samples = {};

            std::ifstream replicate_names{ i_opts.in_replicates_fname };
            std::ofstream averages{ i_opts.out_avgs_fname };

            // Create map of sample names with the base sample name as the key
            // And a list of its associated sample names as the value
            for (std::string current_line; std::getline(replicate_names, current_line); )
                {
                    // Get the base sample name and remove it from the current line
                    std::string delimiter = "\t";
                    size_t pos = current_line.find(delimiter);
                    std::string base_sample = current_line.substr(0, pos);
                    current_line.erase(0, pos + delimiter.length());

                    // Loop through the rest of the sample names on the line
                    // And add them to a vector
                    std::vector<std::string> sample_list;
                    while( ( pos = current_line.find( delimiter ) ) != std::string::npos )
                        {
                            std::string sample = current_line.substr( 0, pos );
                            sample_list.emplace_back( sample );
                            current_line.erase(0, pos + delimiter.length());
                        }
                    // Add the last sample to the vector
                    sample_list.emplace_back( current_line.substr( 0, current_line.find( "\n" ) ) );

                    name_file_samples.emplace( std::make_pair( base_sample, sample_list ) );
                }

            // Check that sample names in names file match those in the input file
            bool invalid_sample_found;
            std::vector<std::string> invalid_samples = {};
            std::vector<std::string> found_samples = {};
            for ( int i = 0; i < scores.sample_names.size(); i++ )
                {
                    invalid_sample_found = true;
                    //std::cout << "SAMPLE TO FIND: \t" << scores.sample_names[i] << std::endl;

                    // Loop through sample names found in name file
                    for ( auto samples : name_file_samples )
                        {
                            //std::cout << samples.first << std::endl;
                            for ( std::string sample : samples.second )
                                {
                                    //std::cout << sample << "\t" << scores.sample_names[i] << "\n";

                                    // If samples in input & name files match, write to output file
                                    if ( std::find( scores.sample_names.begin(),
                                                    scores.sample_names.end(),
                                                    sample ) != scores.sample_names.end() )
                                        {
                                            invalid_sample_found = false;
                                            found_samples.emplace_back( sample );
                                        }
                                    else if ( boost::algorithm::find_backward( invalid_samples.begin(),
                                                                               invalid_samples.end(),
                                                                               samples.first ) == invalid_samples.end() )
                                        {
                                            invalid_samples.emplace_back( samples.first );
                                            std::cout << "Warning: invalid sample name found for: " << samples.first << std::endl;
                                        }
                                }
                        }
                }

            // Write base sample names as row headers in output file
            averages << "Sequence name";
            for(auto sample: name_file_samples)
                {
                    if ( std::find( invalid_samples.begin(), invalid_samples.end(), sample.first ) == invalid_samples.end() )
                        {
                            averages << "\t" << sample.first;
                        }
                }
            averages << "\n";

            // TODO: implent "duplicate sample" warning
            bool duplicate_samples_found = false;
            std::vector<std::string> duplicate_samples;
            int scores_found = 0;
            for ( int pep_index = 0; !invalid_sample_found && pep_index < scores.pep_names.size(); pep_index++ )
                {
                    found_samples = {};
                    averages << scores.pep_names[pep_index];
                    for( int sample_index = 0; sample_index < scores.sample_names.size(); sample_index++)
                        {
                            for( auto sample: name_file_samples )
                                {
                                    if( scores.sample_names[sample_index].find(sample.first) != std::string::npos)
                                        {
                                            // Check that current sample is not a duplicate sample; skip over it if so
                                            if( boost::algorithm::find_backward( found_samples.begin(), found_samples.end(), scores.sample_names[sample_index] )
                                                != found_samples.end() )
                                                {
                                                    // Add sample to list of duplicate samples to be printed in warning;
                                                    // Ensure that it is only added once
                                                    if( boost::algorithm::find_backward( duplicate_samples.begin(), duplicate_samples.end(), scores.sample_names[sample_index])
                                                        == duplicate_samples.end() )
                                                        {
                                                            duplicate_samples_found = true;
                                                            duplicate_samples.emplace_back( scores.sample_names[sample_index] );
                                                        }
                                                    break;
                                                }
                                            else
                                                {
                                                    found_samples.emplace_back( scores.sample_names[sample_index] );
                                                }

                                            sample_map[sample.first].emplace_back( scores.scores.at(sample_index, pep_index) );
                                            scores_found++;
                                            break;
                                        }
                                }
                        }

                    // Check that each sample had a score associated with it; throw error if not
                    if ( scores_found < scores.sample_names.size() )
                        {
                            std::cout << "Error: missing sample score" << std::endl;
                        }

                    float rep_total;
                    float rep_avg = 0.0f;
                    for( auto sample: name_file_samples )
                        {
                            if ( std::find( invalid_samples.begin(),
                                            invalid_samples.end(),
                                            sample.first ) == invalid_samples.end() )
                                  {
                                      // Add all the replicate values for a given sample, then find its average
                                      rep_total = 0;
                                      for( int rep_val: sample_map[sample.first] )
                                          {
                                               rep_total += rep_val;
                                          }
                                      rep_avg = rep_total / (float) sample_map[sample.first].size();
                                      averages << "\t" << rep_avg;
                                  }
                        }
                    averages << "\n";

                    // Reset the sequence map
                    for( auto sample: name_file_samples )
                        {
                            sample_map[sample.first].clear();
                        }
                }

            // Print warning with all duplicate samples found
            if ( duplicate_samples_found )
            {
                std::cout << "Warning: duplicate samples in input file: ";
                for (std::string sample : duplicate_samples)
                {
                    std::cout << sample;
                }
                std::cout << "\n";
            }
        }
}
