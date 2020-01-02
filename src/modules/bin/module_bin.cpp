#include <numeric>

#include "module_bin.h"
#include "omp_opt.h"
#include "peptide_bin.h"
#include "time_keep.h"

module_bin::module_bin() = default;

void module_bin::run( options *opts )
{
    options_bin *b_opts = (options_bin*) opts;
    time_keep::timer timer;
    omp_set_num_threads( 4 );

    timer.start();

    peptide_score_data_sample_major input_data;
    peptide_scoring::parse_peptide_scores( input_data,
                                           b_opts->input_scores_fname
                                         );

    auto peptide_counts = sum_counts( input_data.scores );

    timer.stop();


    std::cout << "Took " << timer.get_elapsed() << " seconds.\n";
}

std::vector<double>
module_bin::sum_counts( const labeled_matrix<double,std::string>& data )
{
    std::vector<double> ret_counts;

    for( std::size_t row_idx = 0; row_idx < data.nrows(); ++row_idx )
        {
            double sum = std::accumulate( data.row_begin( row_idx ),
                                          data.row_end( row_idx ),
                                          0.0
                                        );

            ret_counts.emplace_back( sum );
        }
    return ret_counts;
}

