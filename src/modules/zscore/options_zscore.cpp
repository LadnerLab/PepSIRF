#include "options_zscore.h"
#include <sstream>

options_zscore::options_zscore() = default;

std::string options_zscore::get_arguments()
{
    std::ostringstream str_stream;

    str_stream
        << "--scores      " << in_fname << "\n"
        << "--bins        " << in_bins_fname << "\n"
        << "--trim        " << trim_percent << "\n"
        << "--hdi         " << hpd_percent << "\n"
        << "--output      " << out_fname << "\n"
        << "--nan_report  " << nan_report_fname << "\n"
        << "--num_threads " << num_threads << "\n"
        << "--logfile     " << logfile << "\n\n";

    return str_stream.str();
}

