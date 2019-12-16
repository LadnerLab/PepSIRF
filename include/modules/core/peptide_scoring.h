#ifndef PEPTIDE_SCORING_HH_INCLUDED
#define PEPTIDE_SCORING_HH_INCLUDED

#include <vector>
#include <string>

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

#endif // PEPTIDE_SCORING_HH_INCLUDED
