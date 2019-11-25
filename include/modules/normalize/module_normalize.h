#ifndef MODULE_NORMALIZE_HH_INCLUDED
#define MODULE_NORMALIZE_HH_INCLUDED

#include <string>
#include <vector>

#include "module.h"
#include "options.h"

/**
 * Alias struct peptide_score_data to a type that is 'sample major' (column major with 
 * respect to the input file)
 **/
typedef struct peptide_score_data peptide_score_data_sample_major;

/**
 * Alias struct peptide_score_data to a type that is 'peptide major' (row major with 
 * respect to the input file)
 **/

typedef struct peptide_score_data peptide_score_data_peptide_major;

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
     * using the column_sum method. For this method, 
     * the sum of each column is divided by 1e6.
     * @param cols A vector of doubles, where each value is the 
     *        sum of a column
     * @post Each value in cols is now 1/1e6 its original value
     **/
    void norm_counts_col_sum( std::vector<double>& cols );

};

/**
 * A struct to store peptide score data.
 * This struct is designed to represent 
 * all of the data present in a file that is output by 
 * the demultiplexing module.
 **/
struct peptide_score_data
{
    /**
     * A vector of vectors of scores. 
     * If this struct is peptide-major, scores[ x ][ y ] will 
     * return the score of peptide x in sample y.
     **/
    std::vector<std::vector<double>> scores;

    /**
     * The names of the peptides, in order in which they were
     * found. 
     **/
    std::vector<std::string> pep_names;

    /**
     * The names of the samples, in order in which 
     * they were found.
     **/
    std::vector<std::string> sample_names;
};


#endif // MODULE_NORMALIZE_HH_INCLUDED
