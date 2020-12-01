#ifndef MODULE_P_ENRICH_HH_INCLUDED
#define MODULE_P_ENRICH_HH_INCLUDED
#include <string>
#include <type_traits>

#include "module.h"
#include "options_p_enrich.h"
#include "peptide_scoring.h"
#include "paired_score.h"
#include <map>

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

    /**
     * Determine if the threshold for a given pair of values
     * is met.
     * For a given pair of values to meet a given threshold,
     * the smallest value must be greater than the smallest threshold
     * and the largest value must be greater than the largest threshold.
     * @param values The values to test
     * @param threshold the threshold value that must be met.
     * @returns true if the threshold is met, false otherwise.
     **/
    template<typename ValType>
    bool thresholds_met( const std::pair<ValType,ValType> values,
                             const std::vector<double> thresholds
                           )
    {
        const ValType& a = values.first;
        const ValType& b = values.second;

        if( std::max( a, b ) >= *std::max_element( thresholds.begin(), thresholds.end() )
            && std::min( a, b ) >= *std::min_element( thresholds.begin(), thresholds.end() )
          )
            return true;
        else
            return false;
    }

    /**
     * Get the sums for raw scores in a sequence of
     * paired_scores.
     * @note Statically checks that Iterator iterates over
     *       values of type paired_score.
     * @tparam Iterator an iterator with operator++ over wich
     *         to evaluate items.
     * @param begin The first item in the range [begin, end)
     * @param end the last item in the range [begin, end)
     * @returns a pair of doubles, the first being the sum of each
     *          item's raw_score.first value. The second value
     *          is similarly the sum of each item's raw_score.second
     **/
    template<typename Iterator>
    std::pair<double,double> get_raw_sums( Iterator begin,
                                           Iterator end
                                         )
    {
        std::pair<double,double> sums{ 0.0, 0.0 };

        for( std::map<std::string,paired_score>::iterator current = begin; current != end; ++current )
            {
                sums.first  += current->second.raw_score.first;
                sums.second += current->second.raw_score.second;
            }

        return sums;
    }

    /**
     * Get enrichment candidates for peptides in a sample pair.
     * @param enrichment_candidates pointer to map of paired score values for
     *        each peptide where the key is the petide name
     * @param matrix_score_data pointer to data containing normalized scores.
     * @param raw_score_data pointer to data containing raw scores. If this is
     *        a nullptr, raw scores will be written as zero.
     * @param sample_names the pair of samples to get enrichment candidates for.
     * @param returns void
     **/
    void get_enrichment_candidates( std::map<std::string,paired_score> *enrichment_candidates,
                                    const peptide_score_data_sample_major *matrix_score_data,
                                    const peptide_score_data_sample_major *raw_score_data,
                                    const std::pair<std::string,std::string> sample_names
                                  );
};

#endif // MODULE_P_ENRICH_HH_ENCLUDED
