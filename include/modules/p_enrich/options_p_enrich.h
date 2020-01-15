#ifndef OPTIONS_P_ENRICH_HH_INCLUDED
#define OPTIONS_P_ENRICH_HH_INCLUDED
#include "options.h"

#include <string>

class options_p_enrich : public options
{
 public:

    options_p_enrich()
        : DEFAULT_OUTPUT_DIRNAME{ "paired" },
          DEFAULT_OUT_FNAME_JOIN{ "~" }
          {}

    using score_type = double;
    using score_pair = std::pair<score_type,score_type>;

    std::string get_arguments();

    std::string in_zscore_fname;
    score_pair zscore_params;

    std::string in_norm_scores_fname;
    score_pair norm_scores_params;

    std::string in_raw_scores_fname;
    score_pair raw_scores_params;

    std::string out_dirname;

    std::string out_fname_join;

    std::string out_suffix;

    const std::string DEFAULT_OUTPUT_DIRNAME;
    const std::string DEFAULT_OUT_FNAME_JOIN;

};

#endif // OPTIONS_P_ENRICH_HH_INCLUDED
