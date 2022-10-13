#include <iostream>
#include <fstream>
#include <numeric>
#include <iomanip>

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
            std::unordered_map<std::string, std::vector<int>> sequence_map = {};
        
            std::ifstream replicate_names{ i_opts.in_replicates_fname };
            std::ofstream averages{ i_opts.out_avgs_fname };

            averages << "Sequence name";

            // Get sequence names from the name file
            std::string current_line;
            std::vector<std::string> lines_from_file;
            std::vector<std::string> split_line;
            while( std::getline( replicate_names, current_line ) )
                {
                    lines_from_file.emplace_back(current_line);
                }
        
            // Create list of base sequence names
            std::vector<std::string> base_sequence_names;
            std::string delimiter = "\t";
            std::string sequence;
            std::vector<int> replicate_vector;
            for( std::string line: lines_from_file )
                {
                    // Get base sequence from current line
                    sequence = line.substr(0, line.find(delimiter));

                    // Add to base sequence list, write to file, and place in map as key
                    base_sequence_names.emplace_back(sequence);
                    averages << "\t" << sequence;
                    sequence_map.emplace(std::make_pair(sequence, replicate_vector));
                }

            averages << "\n";

            for( int pep_index = 0; pep_index < scores.pep_names.size(); pep_index++ )
                {
                    averages << scores.pep_names[pep_index];
                    for( int sample_index = 0; sample_index < scores.sample_names.size(); sample_index++)
                        {
                            for( std::string base_seq: base_sequence_names )
                                {
                                    if( scores.sample_names[sample_index].find(base_seq) != std::string::npos)
                                        {
                                            sequence_map[base_seq].emplace_back( scores.scores.at(sample_index, pep_index) );
                                            break;
                                        }
                                }
                        }

                    float rep_total;
                    float rep_avg = 0.0f;
                    for( std::string base_seq: base_sequence_names )
                        {
                            rep_total = 0;
                            for( int rep_val: sequence_map[base_seq] )
                                {
                                     rep_total += rep_val;
                                }
                            rep_avg = rep_total / (float) sequence_map[base_seq].size();
                            averages << "\t" << rep_avg;
                        }
                    averages << "\n";
                    
                    // Reset the sequence map
                    for( std::string base_seq: base_sequence_names )
                        {
                            sequence_map[base_seq].clear();
                        }
                }
            
        }


}
