#include <numeric>

#include "logger.h"
#include "module_bin.h"
#include "omp_opt.h"
#include "time_keep.h"
#include <fstream>
#include <iomanip>

module_bin::module_bin()
{
    name = "Bin";
}

void module_bin::run( options *opts )
{
    options_bin *b_opts = (options_bin*) opts;
    time_keep::timer timer;
    int min_bin_size_extension = 1;
    timer.start();

    // testing purposes -- remove
    info("Bin module has started!");

    peptide_score_data_sample_major input_data;
    peptide_scoring::parse_peptide_scores( input_data,
                                           b_opts->input_scores_fname
                                         );

    auto peptide_counts = sum_counts( input_data.scores.transpose() );
    std::unordered_map<std::string,double> scored_probes;

    // associate each probe with its score
    for( std::vector<double>::size_type idx = 0;
         idx < peptide_counts.size();
         ++idx
       )
        {
            scored_probes.emplace( std::make_pair( input_data.pep_names[ idx ],
                                                   peptide_counts[ idx ]
                                                 )
                                 );
        }

    auto ranked_probes = rank_probes( scored_probes,
                                      b_opts->rounding_factor
                                    );

    auto probe_bins = bin_ranked_probes( ranked_probes,
                                         b_opts->min_bin_size
                                       );
    // while min bin size is greater then smallest bin size, increase min bin size
    while( b_opts->min_bin_size > probe_bins.smallest().size() )
        {
            probe_bins = bin_ranked_probes( ranked_probes,
                            b_opts->min_bin_size + min_bin_size_extension
                                        );
            min_bin_size_extension++;
        }

    std::ofstream output_file( b_opts->output_bins_fname,
                               std::ios_base::out
                             );

    peptide_bin_io::write_bins( output_file, probe_bins );

    timer.stop();

    std::cout << "Took " << timer.get_elapsed() << " seconds.\n";
}

bin_collection module_bin::bin_ranked_probes( const probe_rank& ranked_probes,
                                              const std::size_t min_size
                                            )
{
    std::vector<double> keys;
    keys.reserve( ranked_probes.get_probe_ranks().size() );

    for( const auto& kv_pair : ranked_probes.get_probe_ranks() )
        {
            keys.emplace_back( kv_pair.first );
        }

    std::sort( keys.begin(),
               keys.end(),
               []( const double a, const double b ) -> bool
               { return a > b; }
             );

    bin_collection bins;

    bins.add_bin( peptide_bin() );
    auto current_bin = bins.end() - 1;

    for( const auto key : keys )
        {
            const auto& current_rank_peptides = ranked_probes
                                                .get_probe_ranks()
                                                .find( key )->second;

            current_bin->add_peptides( current_rank_peptides.begin(),
                                       current_rank_peptides.end()
                                     );

            if( current_bin->size()
                >= min_size && key != *( keys.end() - 1 )
              )
                {
                    bins.add_bin( peptide_bin() );
                    current_bin = bins.end() - 1;
                }
        }

    return bins;
}


std::vector<double>
module_bin::sum_counts( const labeled_matrix<double,std::string>& data )
{
    std::vector<double> ret_counts;
    ret_counts.reserve( data.nrows() );

    for( std::size_t row_idx = 0; row_idx < data.nrows(); ++row_idx )
        {
            double sum = std::accumulate( data.row_begin( row_idx ),
                                          data.row_end( row_idx ),
                                          static_cast<double>( 0.0 )
                                        );

            ret_counts.emplace_back( sum );
        }
    return ret_counts;
}

probe_rank module_bin::rank_probes( const std::unordered_map<std::string,double>&
                                    probes_with_scores,
                                    const std::size_t rounding_factor
                                  )
{
    probe_rank ret_val{ rounding_factor };

    for( const auto& scored_probe : probes_with_scores )
        {
            ret_val.rank_probe( scored_probe.second,
                                scored_probe.first
                              );
        }
    return ret_val;
}
