#ifndef MODULE_ENRICH_HH_INCLUDED
#define MODULE_ENRICH_HH_INCLUDED
#include <string>
#include <type_traits>

#include "module.h"
#include "options_enrich.h"
#include "peptide_scoring.h"
#include "paired_score.h"
#include <map>

class module_enrich : public module
{
public:
    using sample_type = std::vector<std::string>;

    /**
     * Default constructor
     */
    module_enrich();

    /**
     * Run the bin module, passing
     * the options specified by a user.
     **/
    void run( options *opts );

    /**
     * Parse a list of sample names from a stream.
     * @param file The stream from which samples should be read.
     * @pre file.good()
     * @returns a vector of sets of strings representing the input
     *          sample name sets.
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
    bool pair_threshold_met(
        const std::pair<ValType,ValType> values,
        const std::pair<ValType,ValType> thresholds
    ) {
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

    /**
     * Determine if the threshold for a given list of values
     * is met.
     * For a given list of values to meet a given threshold,
     * the smallest value must be greater than the smallest threshold
     * and the largest value must be greater than the largest threshold.
     * @param values The values to test
     * @param threshold the threshold value that must be met.
     * @returns true if the threshold is met, false otherwise.
     **/
    template<typename ValType>
    bool thresholds_met(
        const std::vector<ValType> values,
        const std::vector<double> thresholds
    ) {
        if(
            *std::max_element( values.begin(), values.end() )
                >= *std::max_element( thresholds.begin(), thresholds.end() )
            && *std::min_element( values.begin(), values.end() )
                >= *std::min_element( thresholds.begin(), thresholds.end() )
        )
            return true;
        else
            return false;
    }

    /**
     * Get the sums for raw scores in a sequence of
     * raw_scores.
     * @param raw_scores a vector containing vectors of raw scores.
     * @returns a vector of doubles, each value is a sample score being the
     *          sum of all raw scores at each consecutive index sample.
     *          Ex. samples assayed in 2 replicate(duplicate) will result in
     *          a vector of 2 raw sums where a 4 replicate will result in
     *          a vector of 4 raw sums.
     **/
    std::vector<double> get_raw_sums( std::vector<std::vector<double>> raw_scores )
    {

        std::vector<double> sums( raw_scores[0].size(), 0.0 );
        std::size_t raw_score_lists_idx;

        for(
            raw_score_lists_idx = 0;
            raw_score_lists_idx < raw_scores.size();
            raw_score_lists_idx++ 
        ) {
            std::transform(
                raw_scores[raw_score_lists_idx].begin(),
                raw_scores[raw_score_lists_idx].end(),
                sums.begin(),
                sums.begin(),
                std::plus<double>()
            );
        }

        return sums;
    }

    /**
     * Get raw scores for each peptide in specified sample.
     * @param raw_scores_dest The pointer to the location where raw score vectors will be stored.
     * @param raw_score_data The pointer to the raw score data used to find raw scores.
     * @param sample_names The pointer to the sample name vectors to get raw scores from.
     **/
    void get_raw_scores(
        std::vector<std::vector<double>> *raw_scores_dest,
        const peptide_score_data_sample_major *raw_score_data,
        const std::vector<std::string> sample_names
    );

    /**
     * Get enrichment candidates for peptides in a sample list. Each candidate is a pair:
     * (peptide name, list of scores for specified sample columns)
     * @param enrichment_candidates pointer to map of score(s) for
     *        each peptide where the key is the peptide name and the value is a list of >=1 scores.
     * @param matrix_score_data pointer to data containing normalized scores.
     * @param sample_names the list of samples to get enrichment candidates for.
     * @param returns void
     **/
    void get_enrichment_candidates(
        std::map<std::string,std::vector<double>> *enrichment_candidates,
        const peptide_score_data_sample_major *matrix_score_data,
        const std::vector<std::string> sample_names
    );
};


#endif /* MODULE_ENRICH_HH_ENCLUDED */

