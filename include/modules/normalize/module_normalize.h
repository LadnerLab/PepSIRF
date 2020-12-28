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
     * Finding the average for each peptide in control data. Row/peptide major,
     * all samples will be iterated through, but only samples found in neg_filter
     * will be added in calculation of mean for the current peptide.
     * @param control The pointer to data that will be analysed to find scores for
     *        mean calculation.
     * @param neg_filter The pointer to a set containing all sample names used in
     *        calculation.
     * @param pep_averages The pointer to the destination where the averages for
     *        each peptide will be stored (ordered).
     **/
    void get_neg_average( peptide_score_data_sample_major *control,
                                            std::unordered_set<std::string> *neg_filter,
                                            std::vector<double> *pep_averages );


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
     *
     **/
    void filter_neg_control_start( peptide_score_data_sample_major *score_data,
                              std::unordered_set<std::string> *neg_filter,
                              std::string start );

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
     * @param peptide_avgs The means of provided negative controls.
     **/
    void compute_diff( matrix<double> *norm_diffs,
                        std::vector<double> *peptide_avgs );

    /**
     * Compute the difference and ratio:
     * [(individual sample count - negative average) / negative average]
     * for each sample column in a set of data.
     * @param norm_diff_ratios Pointer to the location to store the difference
     *        ratios that were found, one per column (sample).
     * @param peptide_avgs The means of provided negative controls.
     **/
    void compute_diff_ratio( matrix<double> *norm_diff_ratios,
                        std::vector<double> *peptide_avgs );

    /**
     * Compute the ratio of norm socres to negative control score averages for
     * each sample column in a set of data.
     * @param norm_diff_ratios Pointer to the location to store the difference
     *        ratios that were found, one per column (sample).
     * @param peptide_avgs The means of provided negative controls.
     **/
    void compute_ratio( matrix<double> *norm_diff_ratios,
                        std::vector<double> *peptide_avgs );
};

#endif // MODULE_NORMALIZE_HH_INCLUDED
