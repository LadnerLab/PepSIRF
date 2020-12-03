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
    void get_sum( std::vector<double> *dest,
                  matrix<double> *src
                );

    /**
     * Normalize the counts in a vector of column sums
     * by a constant factor.
     * @param cols A vector of doubles, where each value is the
     *        sum of a column
     * @param factor The factor by which to normalize values in cols
     * @post Each value in cols is now 1/factor its original value
     **/
    void constant_factor_normalization( std::vector<double> *cols,
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
    void normalize_counts( matrix<double> *original_scores,
                           const std::vector<double> *norm_factors
                         );

    /**
     * Compute the size factors for each column in a set of data.
     * This method implements equation 5 in Anders and Huber 2010.
     * @param size_factors The location to store the size factors that
     *        were found, one per column (sample).
     * @param data The counts to get the size factors of.
     **/
    void compute_size_factors( std::vector<double> *size_factors,
                               const matrix<double> *data
                             );

    /**
     * Compute the difference between norm scores and negative control score
     * average for each sample column in a set of data.
     * @param norm_diffs Pointer to the location to store the differences that
     *        were found, one per column (sample).
     * @param avg_neg_score The mean of provided negative control score.
     * @param data The counts to get the diffs for.
     **/
    void compute_neg_diff( matrix<double> *norm_diffs,
                           peptide_score_data_sample_major *neg_scores,
                           const matrix<double> *data );


    /**
     * Compute the difference and ratio:
     * [(individual sample count - negative average) / negative average]
     * for each sample column in a set of data.
     * @param norm_diff_ratios Pointer to the location to store the difference
     *        ratios that were found, one per column (sample).
     * @param avg_neg_score The mean of provided negative control score.
     * @param data The counts to get the diff ratios for.
     **/
    void compute_neg_diff_ratio( matrix<double> *norm_diff_ratios,
                                 peptide_score_data_sample_major *neg_scores,
                                 const matrix<double> *data );
};

#endif // MODULE_NORMALIZE_HH_INCLUDED
