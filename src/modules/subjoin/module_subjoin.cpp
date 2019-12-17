#include "module_subjoin.h"
#include "matrix.h"

#include <iostream>
module_subjoin::module_subjoin() = default;

void module_subjoin::parse_namelist( std::vector<std::string>& dest,
                                     std::istream& file
                                   )
{
    std::string line;
    while( std::getline( file, line ).good() )
        {
            dest.push_back( line );
        }
}

void module_subjoin::run( options *opts )
{
    options_subjoin *s_opts = (options_subjoin*) opts;
    std::ifstream names_list( s_opts->names_list_fname,
                              std::ios_base::in
                            );

    std::vector<std::string> peptide_name_list;

    parse_namelist( peptide_name_list, names_list );

    peptide_score_data_sample_major original_scores;
    peptide_scoring::parse_peptide_scores( original_scores,
                                           s_opts->in_matrix_fname
                                         );

    original_scores.scores = original_scores.scores.filter_rows( peptide_name_list );
    original_scores.pep_names = peptide_name_list;
    peptide_scoring::write_peptide_scores( s_opts->out_matrix_fname, original_scores );

}
