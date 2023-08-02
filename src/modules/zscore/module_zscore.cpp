#include "logger.h"
#include "module_zscore.h"
#include "time_keep.h"
#include "stats.h"
#include "omp_opt.h"
#include "noop_iterator.h"
#include "file_io.h"

#include <map>
#include <limits>
#include <iostream>
#include <fstream>


module_zscore::module_zscore()
{
    name = "Zscore";
}

void module_zscore::run( options *opts )
{
    std::map<std::string, std::size_t> names_map;
    std::size_t names_size = 0;
    options_zscore *z_opts = (options_zscore*) opts;
    time_keep::timer time;
    time.start();

    Log::info("Zscore module has started!\n");

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


    for(std::size_t index = 0; index < peptide_bins.size(); ++index)
        {
            for(auto peptidename : peptide_bins.bins[index].peptide_names)
                {
                    ++names_map[peptidename];
                    ++names_size;
                }
        }

    if(names_size != input.pep_names.size())
        {
            throw std::runtime_error("Bins file does not match peptide file");
        }

    for(size_t index = 0; index < input.pep_names.size(); ++index)
        {
            if(names_map[input.pep_names[index]] != 1)
            {
                throw std::runtime_error("Bins file does not match peptide file");
            }
        }


    auto& zscore_matrix = input.scores;
    omp_set_num_threads( z_opts->num_threads );
    std::vector<nan_report> nan_values;

    // for each sample
    for( const auto &sample_pair : zscore_matrix.get_row_labels() )
        {
            std::string sample_name = sample_pair.first;

            calculate_zscores( peptide_bins,
                               z_opts,
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
