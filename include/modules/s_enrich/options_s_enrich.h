#ifndef OPTIONS_S_ENRICH_HH_INCLUDED
#define OPTIONS_S_ENRICH_HH_INCLUDED
#include "options.h"
#include <string>

class options_s_enrich : public options
{
 public:

    options_s_enrich();

    std::string get_arguments();

    /**
     * Name of the directory to write 
     * output files to.
     **/
    std::string out_dirname;

    /**
     * Suffix to append to output files.
     **/
    std::string out_suffix;

    /**
     * Name of the score matrix containing zscores 
     * for peptides in samples
     **/
    std::string in_zscore_fname;

    /**
     * Name of the file containing normalized 
     * scores for peptides.
     **/
    std::string in_norm_score_fname;

    /**
     * Name of the file containing raw count data.
     **/
    std::string in_raw_count_fname;

    /**
     * The minimum zscore a peptide must 
     * have in a sample to be considered 
     * enriched
     **/
    double min_zscore;

    /**
     * The minimum normalized score a peptide must 
     * have in a sample to be considered 
     * enriched
     **/
    double min_norm_score;

    /**
     * The minimum raw count a peptide must 
     * have in a sample to be considered 
     * enriched
     **/
    double min_raw_count;
};


#endif // OPTIONS_S_ENRICH_HH_INCLUDED
