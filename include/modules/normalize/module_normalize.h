#ifndef MODULE_NORMALIZE_HH_INCLUDED
#define MODULE_NORMALIZE_HH_INCLUDED

#include <string>
#include <vector>
#include <cmath>

#include "module.h"
#include "options.h"
#include "peptide_scoring.h"

/**
 * Class for the PepSIRF normalization module.
 **/
class module_normalize : public module
{

 public:
    std::string name;

    module_normalize();

    std::string get_name();

    void run( options *opts );

    /**
     * Get the sum of each column in a matrix src.
     * @param dest The location to store the sums.
     *        After this function, the n'th position in dest will 
     *        hold the sum of the n'th column of src.
     * @param N x M matrix representing the values parsed from 
     *        a score matrix file.
     **/
    void get_sum( std::vector<double>& dest,
                  matrix<double>& src
                );


    /**
     * Normalize the counts in a vector of column sums 
     * by a constant factor.
     * @param cols A vector of doubles, where each value is the 
     *        sum of a column
     * @param factor The factor by which to normalize values in cols
     * @post Each value in cols is now 1/factor its original value
     **/
    void constant_factor_normalization( std::vector<double>& cols,
                                        const std::size_t factor
                                      );

    /**
     * Normalize the counts in 'original_scores'.
     * Each score in a column will be divided by the respective index in 
     * norm_factors, i.e. each value in column n will be divided by 
     * norm_factors[ n ]. 
     * @param original_scores A score matrix of the original scores.
     * @param norm_factors The column normalization factors to adjust the 
     *        score of each in original_score by.
     **/
    void normalize_counts( matrix<double>&
                           original_scores,
                           const std::vector<double>& norm_factors
                         );

    /**
     * Compute the geometric mean of a set of data of type T.
     * @pre T must be a numeric type (std::size_t, double, int, etc.)
     * @pre all of the values in data are greater than zero
     * @note Only non-zero values are considered in the calculation, so 
     *       the 'size' of data is considered the number of non-zero values in 
     *       data.
     * @param data The dataset whose geometric mean to compute.
     **/
    template <typename Iterator
             >
        double geom_mean( const Iterator& begin,
                          const Iterator& end
                        )
        {
            std::size_t n = 0;
            double log_sum = 0;

            for( auto index = begin; index != end; ++index )
                {
                    log_sum += std::log( *index == 0 ? 1 : *index );
                    n += *index != 0;
                }

            log_sum /= (double) n;

            return n == 0 ? 0 : std::exp( log_sum );
        }

    /**
     * Compute the size factors for each column in a set of data.
     * This method implements equation 5 in Anders and Huber 2010.
     * @param size_factors The location to store the size factors that 
     *        were found, one per column (sample).
     * @param data The counts to get the size factors of.
     **/
    void compute_size_factors( std::vector<double>& size_factors,
                               const matrix<double>& data
                             );


};

#endif // MODULE_NORMALIZE_HH_INCLUDED
