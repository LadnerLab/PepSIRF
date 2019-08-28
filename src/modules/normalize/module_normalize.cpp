#include <iostream>
#include <omp.h>
#include <algorithm>
#include <fstream>

#include "module_normalize.h"
#include "options_normalize.h"

module_normalize::module_normalize()
{
    name = "Normalize";
}

std::string module_normalize::get_name()
{
    return name;
}

void module_normalize::run( options *opts )
{
    options_normalize *n_opts = (options_normalize*) opts;

    std::string scores_fname = n_opts->peptide_scores_fname;

    omp_set_num_threads( n_opts->num_threads );

    peptide_score_data_sample_major original_scores;

}

void module_normalize::parse_peptide_scores( peptide_score_data_sample_major& dest,
                                             std::string ifname
                                           )
{

}

