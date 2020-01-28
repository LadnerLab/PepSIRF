#ifndef OPTIONS_ZSCORE_HH_INCLUDED
#define OPTIONS_ZSCORE_HH_INCLUDED
#include <string>
#include "options.h"

class options_zscore : public options
{

 public:
    options_zscore();
    std::string get_arguments();

    std::string in_fname;
    std::string in_bins_fname;
    std::string out_fname;
    std::string nan_report_fname;

    double trim_percent;
};

#endif // OPTIONS_ZSCORE_HH_INCLUDED
