#ifndef MODULE_BIN_HH_INCLUDED
#define MODULE_BIN_HH_INCLUDED
#include <vector>
#include <unordered_map>
#include <string>

#include "peptide_scoring.h"
#include "module.h"
#include "options_bin.h"
#include "probe_rank.h"
#include "peptide_bin.h"

class module_bin : public module
{
 public:
    module_bin();

    /**
     * Run the bin module, passing 
     * the options specified by a user.
     **/
    void run( options *opts );

    /**
     * Sum the counts of the rows in a labeled_matrix.
     * @param data The labeled_matrix whose row values should be summed
     * @returns a vector containing the sum of the counts of the rows in the matrix. 
     *          The n'th entry in this vector is the sum of the n'th row in the data 
     *          matrix.
     **/
    std::vector<double>
    sum_counts( const labeled_matrix<double,std::string>& data );

    /**
     * Rank probes based upon their scores and the rounding factor.
     * @param probes_with_scores An unordered map containing 
     *        probe : score pairs.
     * @param rounding_factor the rounding factor to use when calculating 
     *        the rank of a probe from its score.
     * @returns A probe_rank containing the ranked probes.
     **/
    probe_rank rank_probes( const std::unordered_map<std::string,double>&
                            probes_with_scores,
                            const std::size_t rounding_factor
                          );

    bin_collection bin_ranked_probes( const probe_rank& ranked_probes,
                                      const std::size_t min_size
                                    );

};

#endif // MOUDLE_BIN_HH_INCLUDED
