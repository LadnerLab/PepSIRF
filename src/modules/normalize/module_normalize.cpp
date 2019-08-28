#include <iostream>
#include <omp.h>

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

}
