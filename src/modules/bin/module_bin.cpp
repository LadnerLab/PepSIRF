#include "module_bin.h"
#include "peptide_scoring.h"
#include "omp_opt.h"
#include "peptide_bin.h"
#include "time_keep.h"

module_bin::module_bin() = default;

void module_bin::run( options *opts )
{
    options_bin *b_opts = (options_bin*) opts;
    time_keep::timer timer;

    timer.start();

    peptide_score_data_sample_major input_data;
    peptide_scoring::parse_peptide_scores( input_data,
                                           b_opts->input_scores_fname
                                         );

    input_data.scores = input_data.scores.transpose();


    timer.stop();


    std::cout << "Took " << timer.get_elapsed() << " seconds.\n";
}
