#ifndef MODULE_P_ENRICH_HH_INCLUDED
#define MODULE_P_ENRICH_HH_INCLUDED
#include "module.h"
#include "options_p_enrich.h"
#include "peptide_scoring.h"
#include <string>
#include "paired_score.h"

class module_p_enrich : public module
{

public:

    using sample_type = std::pair<std::string,std::string>;
    /**
     * Run the bin module, passing 
     * the options specified by a user.
     **/
    void run( options *opts );

    /**
     * Parse a list of sample names from a stream. 
     * @param file The stream from which samples should be read.
     * @throws std::runtime_exception if the file is not formatted correctly.
     * @pre file.good()
     * @returns a vector of pairs of strings representing the input 
     *          sample name pairs.
     **/
    std::vector<sample_type> parse_samples( std::istream& file );

    /**
     * Determine if the threshold for a given pair of values 
     * is met. 
     * For a given pair of values to meet a given pair of thresholds,
     * each of the items in 'values' must be at least one of the items in
     * thresholds, independent of order.
     * @param values The values to test
     * @param thresholds the thresholds values must meet.
     * @returns true if the threshold is met, false otherwise.
     **/
    template<typename ValType>
    bool pair_threshold_met( const std::pair<ValType,ValType> values,
                             const std::pair<ValType,ValType> thresholds
                           )
    {
        const ValType& a = values.first;
        const ValType& b = values.second;
        const ValType& c = thresholds.first;
        const ValType& d = thresholds.second;

        // both a and b are at least either of c and
        // not a and b are both less than either one 
        bool at_least = ( a >= c || a >= d ) && ( b >= c || b >= d );
        bool not_less = ( a < c || a < d ) && ( b < c || b < d );


        return at_least && !not_less;
    }

    std::vector<std::string>
    get_enrichment_candidates( const peptide_score_data_sample_major *zscore_data,
                               const peptide_score_data_sample_major *norm_score_data,
                               const peptide_score_data_sample_major *raw_score_data,
                               const std::pair<std::string,std::string> sample_names
                             );

};

#endif // MODULE_P_ENRICH_HH_ENCLUDED
