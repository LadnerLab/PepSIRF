#include "options_s_enrich.h"
#include <sstream>

options_s_enrich::options_s_enrich() = default;


std::string options_s_enrich::get_arguments()
{
    std::ostringstream str_stream;

    str_stream << " --zscores        " << in_zscore_fname << "\n " 
               << " --min_zscore     " << min_zscore << "\n "
               << " --norm_scores    " << in_norm_score_fname << "\n "
               << " --min_norm_score " << min_norm_score << "\n ";

    if( !in_raw_count_fname.empty() )
        {
            str_stream << " --raw_counts     " << in_raw_count_fname << "\n " 
                       << " --min_raw_count  " << min_raw_count << "\n "
                       ;
        }

    str_stream << "\n";

    return str_stream.str();
}
