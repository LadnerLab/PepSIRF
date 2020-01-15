#include "options_p_enrich.h"
#include <sstream>
#include "stream_tools.h"

std::string options_p_enrich::get_arguments()
{
    // enable ADL
    using namespace std;
    
    std::ostringstream stream;
    stream << "--zscores               " << in_zscore_fname << "\n "
           << "--zscore_constraint     " << zscore_params << "\n "
           << "--norm_scores           " << in_norm_scores_fname << "\n "
           << "--norm_score_constraint " << norm_scores_params << "\n "
        ;

    if( !in_raw_scores_fname.empty() )
        {
           stream << "--raw_scores            " << in_raw_scores_fname << "\n "
                  << "--raw_score_constraint  " << raw_scores_params << "\n ";
        }

    stream << "--outfile_suffix        " << out_suffix << "\n "
           << "--output                " << out_dirname << "\n "
           << "\n"
        ;

    return stream.str();
}
