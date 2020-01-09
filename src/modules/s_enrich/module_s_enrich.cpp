#include <iostream>
#include <unordered_set>
#include <stdexcept>

#include "module_s_enrich.h"
#include "time_keep.h"
#include "omp_opt.h"
#include "fs_tools.h"

void module_s_enrich::run( options *opts )
{
    options_s_enrich *e_opts = (options_s_enrich*) opts;
    time_keep::timer timer;
    timer.start();


    peptide_score_data_sample_major zscores;
    peptide_score_data_sample_major norm_scores;
    peptide_score_data_sample_major raw_scores;

    peptide_score_data_sample_major *zscores_ptr = &zscores;
    peptide_score_data_sample_major *norm_scores_ptr = &zscores;
    peptide_score_data_sample_major *raw_scores_ptr = nullptr;

    peptide_scoring::parse_peptide_scores( zscores, e_opts->in_zscore_fname );
    peptide_scoring::parse_peptide_scores( norm_scores, e_opts->in_norm_score_fname );

    bool raw_counts_included = !e_opts->in_raw_count_fname.empty();

    if( raw_counts_included )
        {
            peptide_scoring::parse_peptide_scores( raw_scores,
                                                   e_opts->in_raw_count_fname
                                                 );
            raw_scores_ptr = &raw_scores;
        }

    auto output_path = fs_tools::path( e_opts->out_dirname );
    bool dir_exists = !fs_tools::create_directories( output_path );

    if( dir_exists )
        {
            std::cout << "WARNING: '" << e_opts->out_dirname << "' Exists, any files with " 
                      << "colliding filenames will be overwritten!\n";
        }

    
    timer.stop();
    std::cout << "Took " << timer.get_elapsed() << " seconds.\n";
}

std::vector<peptide_score<std::string>>
module_s_enrich::get_enrichment_candidates( const peptide_score_data_sample_major *zscore_data,
                                            const peptide_score_data_sample_major *norm_score_data,
                                            const peptide_score_data_sample_major *raw_score_data,
                                            const std::size_t sample_idx
                                          )
{
    std::vector<peptide_score<std::string>> ret_val;
    std::size_t pep_idx = 0;

    for( pep_idx = 0; pep_idx < zscore_data->pep_names.size(); ++pep_idx )
        {
            std::string pep_name = zscore_data->pep_names[ pep_idx ];

            double pep_zscore = zscore_data->scores( sample_idx, pep_idx );
            double pep_norm_score = norm_score_data->scores( sample_idx, pep_idx );
            double pep_raw_score = 0;

            if( raw_score_data != nullptr )
                {
                    pep_raw_score = raw_score_data->scores( sample_idx, pep_idx );
                }

            ret_val.emplace_back( peptide_score<std::string>( pep_name,
                                                              pep_zscore,
                                                              pep_norm_score,
                                                              pep_raw_score
                                                            )
                                );
        }

    return ret_val;
}

void module_s_enrich::write_probe_names( std::ostream &stream,
                                         const std::vector<peptide_score<std::string>>& probes
                                       )
{
    for( auto probe = probes.begin();
         probe != probes.end();
         ++probe
       )
        {
            stream << probe->peptide;
            stream << "\n";
        }
}
