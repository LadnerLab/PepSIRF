#ifndef PEPTIDE_SCORING_HH_INCLUDED
#define PEPTIDE_SCORING_HH_INCLUDED

#include <vector>
#include <string>
#include "matrix.h"

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
 * A struct to store peptide score data.
 * This struct is designed to represent 
 * all of the data present in a file that is output by 
 * the demultiplexing module.
 **/
struct peptide_score_data
{
    /**
     * A matrix of scores.
     * If this struct is peptide-major, scores[ x ][ y ] will 
     * return the score of peptide x in sample y.
     **/
    labeled_matrix<double,std::string> scores;

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

    /**
     * The name of the file this peptide_score_data is from.
     **/
    std::string file_name;

    peptide_score_data( const peptide_score_data& other )
    : scores{ other.scores }, pep_names{ other.pep_names },
      sample_names{ other.sample_names}, file_name{ other.file_name }
       {};

    peptide_score_data() = default;
};

namespace peptide_scoring
{

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
     * Write peptide scores from data to an output stream.. 
     * Output will be written in a score matrix where an entry 
     * (x, y) in the score matrix is the the score of peptide x 
     * in sample y.
     * @param output The output stream to write peptide scores to.
     * @param data the peptide score data to write output to
     **/
    void write_peptide_scores( std::ostream& output,
                               peptide_score_data_sample_major& data
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

}; // namespace peptide_scoring

#endif // PEPTIDE_SCORING_HH_INCLUDED
