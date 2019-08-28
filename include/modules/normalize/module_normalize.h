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
    std::vector<std::vector<std::size_t>> scores;

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
