#ifndef OPTIONS_S_ENRICH_HH_INCLUDED
#define OPTIONS_S_ENRICH_HH_INCLUDED
#include "options.h"
#include <string>

class options_s_enrich : public options
{
 public:

    options_s_enrich();

    std::string get_arguments();

    std::string out_dirname;

    std::string in_zscore_fname;

    std::string in_norm_score_fname;

    std::string in_raw_count_fname;

    double min_zscore;
    double min_norm_score;
    double min_raw_count;

};


#endif // OPTIONS_S_ENRICH_HH_INCLUDED
