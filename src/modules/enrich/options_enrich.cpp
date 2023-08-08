#include "options_enrich.h"
#include <sstream>
#include <chrono>
#include <ctime>
#include "stream_tools.h"

std::string options_enrich::get_arguments()
{
    // enable ADL
    using namespace std;

    std::ostringstream stream;
    stream
        << "--samples                    " << in_samples_fname << "\n"
        << "--threshhold_file            " << threshold_fname << "\n";

    if (!in_raw_scores_fname.empty())
    {
        stream
            << "--raw_scores                 " << in_raw_scores_fname << "\n"
            << "--raw_score_constraint       " << raw_scores_params_str << "\n";
    }

    stream
        << "--enrichment_failure_reason  " << out_enrichment_failure << "\n"
        << "--outfile_suffix             " << out_suffix << "\n"
        << "--join_on                    " << out_fname_join << "\n"
        << "--output_filename_truncate   " << "\n"
        << "--output                     " << out_dirname << "\n"
        << "--logfile                    " << logfile << "\n\n";

    return stream.str();
}

