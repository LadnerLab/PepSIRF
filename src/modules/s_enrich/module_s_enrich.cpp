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

    peptide_scoring::parse_peptide_scores( zscores, e_opts->in_zscore_fname );
    peptide_scoring::parse_peptide_scores( zscores, e_opts->in_norm_score_fname );

    timer.stop();

    
    std::cout << "Took " << timer.get_elapsed() << " seconds.\n";
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
