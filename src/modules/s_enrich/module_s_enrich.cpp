#include <iostream>

#include "module_s_enrich.h"
#include "peptide_scoring.h"
#include "time_keep.h"
#include "omp_opt.h"

void module_s_enrich::run( options *opts )
{
    options_s_enrich *e_opts = (options_s_enrich*) opts;
    time_keep::timer timer;
    timer.start();


    peptide_score_data_sample_major zscores;
    peptide_score_data_sample_major norm_scores;

    peptide_scoring::parse_peptide_scores( zscores, e_opts->in_zscore_fname );
    peptide_scoring::parse_peptide_scores( zscores, e_opts->in_norm_score_fname );

    timer.stop();

    
    std::cout << "Took " << timer.get_elapsed() << " seconds.\n";
}
