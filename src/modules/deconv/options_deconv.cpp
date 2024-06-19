#include <algorithm>
#include <sstream>

#include "options_deconv.h"

std::string options_deconv::get_arguments()
{
    std::ostringstream str_stream;

    // returns the 'string' representation of a bool
    auto bool_str = []( bool opt ) -> std::string
        { return opt == true ? "true" : "false"; };

    str_stream
        << "--linked                  " << linked_fname << "\n"
        << "--threshold               " << threshold << "\n"
        << "--enriched                " << enriched_fname << "\n"
        << "--single_threaded         " << bool_str(single_threaded) << "\n"
        << "--scoring_strategy        " << scoring_strategy << "\n"
        << "--score_filtering         " << bool_str(score_filtering) << "\n"
        << "--id_name_map             " << id_name_map_fname << "\n"
        << "--custom_id_name_map_info " << tuple_to_string(custom_id_name_map_info) << "\n"
        << "--score_tie_threshold     " << score_tie_threshold << "\n"
        << "--score_overlap_threshold " << score_overlap_threshold << "\n"
        << "--output                  " << output_fname << "\n"
        << "--outfile_suffix          " << outfile_suffix << "\n"
        << "--scores_per_round        " << orig_scores_dname << "\n"
        << "--peptide_assignment_map  " << species_peptides_out << "\n"
        << "--mapfile_suffix          " << map_suffix << "\n"
        << "--enriched_file_ending    " << enriched_file_ending << "\n"
        << "--logfile                 " << logfile << "\n\n";

    return str_stream.str();
}

//template<typename T>
void options_deconv::set_info(
    std::tuple<std::string, std::string, std::string>& tup,
    std::string info
) {
    const char *COMMA = ",";
    const int  NUM_REQUIRED_ITEMS = 3;

    std::vector<std::string> matches;
    boost::split(matches, info, boost::is_any_of(COMMA));

    if (matches.size() != NUM_REQUIRED_ITEMS)
    {
        Log::error(
            "Incorrect number of comma-separated values provided."
            " Please try ./pepsirf demux -h for help"
        );
    }

    std::get<0>(tup) = matches[0];
    std::get<1>(tup) = matches[1];
    std::get<2>(tup) = matches[2];
}

std::string options_deconv::tuple_to_string(std::tuple<std::string, std::string, std::string>& tup)
{
    std::ostringstream str;

    str << std::get<0>(tup) << "," << std::get<1>(tup) << "," << std::get<2>(tup);

    return str.str();
}
