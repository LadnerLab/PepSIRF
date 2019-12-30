#include "module_zscore.h"
#include "time_keep.h"
#include "stats.h"
#include "omp_opt.h"
#include <iostream>
#include <fstream>


module_zscore::module_zscore() = default;

void module_zscore::run( options *opts )
{
    options_zscore *z_opts = (options_zscore*) opts;
    time_keep::timer time;
    time.start();

    // parse the input
    peptide_score_data_sample_major input;

    peptide_scoring::parse_peptide_scores( input, z_opts->in_fname );

    // transpose the input matrix
    input.scores = input.scores.transpose();

    // calculate the zscores
    for( std::uint32_t row_idx = 0; row_idx < input.scores.nrows(); ++row_idx )
        {
            stats::zscore( input.scores.row_begin( row_idx ),
                           input.scores.row_end( row_idx ),
                           input.scores.row_begin( row_idx )
                         );
        }

    // transpose the output matrix
    input.scores = input.scores.transpose();


    // write output
    peptide_scoring::write_peptide_scores( z_opts->out_fname, input );

    time.stop();

    std::cout << "Took " << time.get_elapsed() << " seconds.\n";
}
void
module_zscore::calculate_zscores( const bin_collection& peptide_bins,
                                  const double trim_percent,
                                  labeled_matrix<double,std::string>& scores,
                                  const std::uint32_t sample_num
                                )
{
    std::vector<double> current_row_counts;

    for( const auto& bin : peptide_bins )
        {
            // get the counts for the peptides in this group
            for( const auto& peptide : bin )
                {
                    current_row_counts.emplace_back( scores( sample_num, peptide ) );
                }

            std::sort( current_row_counts.begin(),
                       current_row_counts.end()
                     );

            std::size_t trim_length = current_row_counts.size() * ( 0.01 * trim_percent );
            auto new_begin = current_row_counts.begin() + trim_length;
            auto new_end = current_row_counts.end() - trim_length;

            double mean = stats::arith_mean( new_begin,
                                             new_end
                                           );

            double stdev = stats::stdev( new_begin,
                                         new_end,
                                         mean
                                       );

            for( const auto& peptide : bin )
                {

                    auto &current_value = scores( sample_num, peptide );
                    double zscore = stats::zscore( current_value, mean, stdev );
                    current_value = zscore;
                }

            current_row_counts.clear();
        }
}
