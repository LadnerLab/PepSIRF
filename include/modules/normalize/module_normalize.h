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
     * Parse the scores of peptides from a file. 
     * Data (scores, peptide names, and sample names) will all 
     * be stored in the appropriate members of dest.
     * @param dest The destination for scores. 
     * @param ifname The name of the input file to parse. This file should be 
     *        in the format output by the 'demux' module. 
     * @note to reduce resource usage and improve efficiency of 
     *       iteration items are stored in dest in sample-major order,
     *       so accessing dest.scores[ x ][ y ] returns the score for 
     *       the y'th peptide in sample x.
     **/
    void parse_peptide_scores( peptide_score_data_sample_major& dest,
                               std::string ifname
                             );

    /**
     * Write peptide scores from data to a file. 
     * Output will be written in a score matrix where an entry 
     * (x, y) in the score matrix is the the score of peptide x 
     * in sample y.
     * @param dest_fname The name of file to write output to
     * @param data the peptide score data to write output to
     **/
    void write_peptide_scores( std::string dest_fname,
                               peptide_score_data_sample_major& data
                             );

    /**
     * Get the sum of each column in a matrix src.
     * @param dest The location to store the sums.
     *        After this function, the n'th position in dest will 
     *        hold the sum of the n'th column of src.
     * @param N x M matrix representing the values parsed from 
     *        a score matrix file.
     **/
    void get_sum( std::vector<double>& dest,
                  std::vector<std::vector<double>>& src
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
    void normalize_counts( std::vector<std::vector<double>>&
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
    template <typename T>
        double geom_mean( const std::vector<T>& data )
        {
            std::size_t n = 0;
            std::size_t index = 0;
            double log_sum = 0;

            for( index = 0; index < data.size(); ++index )
                {
                    log_sum += std::log( std::max( (T) 1, (T)data[ index ] ) );
                    n += data[ index ] != 0;
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
                               const std::vector<std::vector<double>>& data
                             );


};

#endif // MODULE_NORMALIZE_HH_INCLUDED
