#ifndef OPTIONS_SUBJOIN_HH_INCLUDED
#define OPTIONS_SUBJOIN_HH_INCLUDED
#include <vector>
#include <string>
#include <unordered_map>

#include "options.h"

namespace evaluation_strategy
{
    enum duplicate_resolution_strategy
    {
        COMBINE,
        INCLUDE,
        IGNORE
    };

    extern std::unordered_map<std::string,duplicate_resolution_strategy> drs_string_map;
    extern std::unordered_map<duplicate_resolution_strategy,std::string> string_drs_map;

    bool is_valid( const std::string& check );

    duplicate_resolution_strategy from_string( const std::string& str );
    std::string to_string( duplicate_resolution_strategy strategy );
};

/**
 * Options for the subjoin module of the
 * PepSIRF package.
 **/
class options_subjoin : public options
{
 public:

    options_subjoin();
    /**
     * Returns a string of the arguments provided
     * to the module.
     **/
    std::string get_arguments();

    std::vector<std::pair<std::string,std::string>> multi_matrix_name_pairs;
    std::vector<std::pair<std::string,std::string>> input_matrix_name_pairs;
    /**
     * The name of the file to write output to.
     **/
    std::string out_matrix_fname;

    /**
     * Boolean option to determine whether
     * sample names or peptide names should be
     * used.
     **/
    bool use_sample_names;
    bool exclude_names;

    evaluation_strategy::duplicate_resolution_strategy duplicate_resolution_strategy;
};

#endif // OPTIONS_SUBJOIN_HH_INCLUDED
