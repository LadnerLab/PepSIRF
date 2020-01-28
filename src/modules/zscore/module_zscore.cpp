#include "module_zscore.h"
#include "time_keep.h"
#include "stats.h"
#include "omp_opt.h"
#include "noop_iterator.h"
#include "file_io.h"

#include <limits>
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

    std::ifstream bins_file( z_opts->in_bins_fname, std::ios_base::in );

    if( bins_file.fail() )
        {
            throw std::runtime_error( "Unable to open the "
                                      "specified bins file for reading"
                                    );
        }
        
    bin_collection peptide_bins = peptide_bin_io::parse_bins( bins_file );

    auto& zscore_matrix = input.scores;
    omp_set_num_threads( z_opts->num_threads );
    std::vector<nan_report> nan_values;

    // for each sample
    for( const auto sample_pair : zscore_matrix.get_row_labels() )
        {
            std::string sample_name = sample_pair.first;

            calculate_zscores( peptide_bins,
                               z_opts->trim_percent,
                               zscore_matrix,
                               sample_name,
                               std::back_inserter( nan_values )
                             );
        }

    // write output
    peptide_scoring::write_peptide_scores( z_opts->out_fname, input );

    if( !z_opts->nan_report_fname.empty() )
        {
            std::ofstream report_file{ z_opts->nan_report_fname };
            report_file << "Probe name\tSample Name\tBin Number\n";
            pepsirf_io::write_file( report_file,
                                    nan_values.begin(),
                                    nan_values.end(),
                                    "\n"
                                  );
        }

    time.stop();

    std::cout << "Took " << time.get_elapsed() << " seconds.\n";
}
