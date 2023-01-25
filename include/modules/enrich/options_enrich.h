#ifndef OPTIONS_ENRICH_HH_INCLUDED
#define OPTIONS_ENRICH_HH_INCLUDED
#include "options.h"

#include <string>
#include <vector>
#include <set>

class options_enrich : public options
{
 public:

    options_enrich()
        : DEFAULT_OUTPUT_DIRNAME{ "paired" },
          DEFAULT_OUT_FNAME_JOIN{ "~" }
          {}

    using score_type = double;
    using score_pair = std::pair<score_type,score_type>;

    std::string get_arguments();

    std::string in_samples_fname;

    std::string in_raw_scores_fname;
    std::string raw_scores_params_str;
    std::vector<double> raw_scores_params;

    std::string threshold_fname;
    std::vector<std::pair<std::string,std::string>> matrix_thresh_fname_pairs;

    std::string out_dirname;

    std::string out_enrichment_failure;

    std::string out_fname_join;

    std::string out_suffix;

    bool truncate_names;

    bool low_raw_reads;

    const std::string DEFAULT_OUTPUT_DIRNAME;
    const std::string DEFAULT_OUT_FNAME_JOIN;

};

#endif // OPTIONS_ENRICH_HH_INCLUDED
