#include <algorithm>
#include <sstream>

#include "options_deconv.h"

std::string options_deconv::get_arguments()
{
    if (logfile.empty())
    {
        logfile = set_default_log();
    }

    std::ostringstream str_stream;

    // returns the 'string' representation of a bool
    auto bool_str = []( bool opt ) -> std::string
        { return opt == true ? "true" : "false"; };

    str_stream << "--linked                  " << linked_fname  << "\n" <<  // not sure why but this one needs one space fewer
        " --threshold               " << threshold <<  "\n" <<
        " --enriched                " << enriched_fname <<  "\n" <<
        " --single_threaded         " << bool_str( single_threaded ) <<  "\n" <<
        " --scoring_strategy        " << scoring_strategy <<  "\n" <<
        " --score_filtering         " << bool_str( score_filtering ) <<  "\n" <<
        " --id_name_map             " << id_name_map_fname << "\n" <<
        " --score_tie_threshold     " << score_tie_threshold << "\n" <<
        " --score_overlap_threshold " << score_overlap_threshold << "\n" <<
        " --output                  " << output_fname <<  "\n"
        " --outfile_suffix          " << outfile_suffix <<  "\n"
        " --scores_per_round        " << orig_scores_dname <<  "\n"
        " --peptide_assignment_map  " << species_peptides_out <<  "\n"
        " --mapfile_suffix          " << map_suffix <<  "\n"
        " --enriched_file_ending    " << enriched_file_ending << "\n" <<
        " --logfile                 " << logfile << "\n"
        ;

    str_stream << "\n";


    return str_stream.str();
}
