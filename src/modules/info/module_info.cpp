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
                                  << sum
                                  << "\n";
                }

        }


}
