#include "module_p_enrich.h"
#include <numeric>
#include "time_keep.h"
#include "omp_opt.h"
#include "fs_tools.h"
#include "predicate.h"


void module_p_enrich::run( options *opts )
{
    options_p_enrich *p_opts = (options_p_enrich*) opts;
    time_keep::timer timer;
    timer.start();


    peptide_score_data_sample_major zscores;
    peptide_score_data_sample_major norm_scores;
    peptide_score_data_sample_major raw_scores;

    peptide_score_data_sample_major *zscores_ptr = &zscores;
    peptide_score_data_sample_major *norm_scores_ptr = &zscores;
    peptide_score_data_sample_major *raw_scores_ptr = nullptr;

    peptide_scoring::parse_peptide_scores( zscores, p_opts->in_zscore_fname );
    peptide_scoring::parse_peptide_scores( norm_scores, p_opts->in_norm_scores_fname );

    bool raw_counts_included = !p_opts->in_raw_scores_fname.empty();

    auto output_path = fs_tools::path( p_opts->out_dirname );
    bool dir_exists = !fs_tools::create_directories( output_path );

    if( dir_exists )
        {
            std::cout << "WARNING: the directory '" << p_opts->out_dirname
                      << "' exists, any files with " 
                      << "colliding filenames will be overwritten!\n";
        }

    using sample_name_set = std::unordered_set<std::string>;
    sample_name_set zscore_sample_names{ zscores.sample_names.begin(),
            zscores.sample_names.end()
            };

    sample_name_set norm_score_sample_names{ norm_scores.sample_names.begin(),
            norm_scores.sample_names.end()
            };

    sample_name_set raw_score_sample_names{ raw_scores.sample_names.begin(),
            raw_scores.sample_names.end()
            };

    // ensure zscore_sample_names and norm_score_sample_names are
    // equal. If included, raw_count_sample_names must also equal
    // zscore equal names, otherwise we do not care.
    if( zscore_sample_names != norm_score_sample_names
        && ( !raw_counts_included
             || ( zscore_sample_names == raw_score_sample_names
                  )
             )
        )
        {
            throw std::runtime_error( "The samplenames provided in each input file are "
                                      "not the same"
                                    );
        }


}
