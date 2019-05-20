#ifndef OPTIONS_HH_INCLUDED
#define OPTIONS_HH_INCLUDED
#include <string>

/*! Data class to contain and handle 
 * arguments that will be passed in from the 
 * command-line.
*/
class options
{
public:
    options(): DEFAULT_NUM_THREADS( 2 ){} //!< Default constructor.

    std::string input_r1_fname; //!< Filename for forward reads, can be a .zip archive or a regular fastq file.
    std::string input_r2_fname; //!< Filename for reverse reads, can be a .zip archive or regular fastq file.
    std::string library_fname; //!< Filename containing a FASTA file containing a library of amino acid peptide sequences.

    int num_threads; //!< The number of threads to use for computation.

    const int DEFAULT_NUM_THREADS; //!< The default number of threads to use
};

#endif //OPTIONS_HH_INCLUDED
