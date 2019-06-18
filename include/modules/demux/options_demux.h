#ifndef OPTIONS_DEMUX_HH_INCLUDED
#define OPTIONS_DEMUX_HH_INCLUDED
#include <string>
#include <sstream>
#include "options.h"

/*! Data class to contain and handle 
 * arguments that will be passed in from the 
 * command-line.
*/
class options_demux: public options
{
public:
    options_demux(): DEFAULT_READ_PER_LOOP( 100000 ), DEFAULT_OUTPUT_FNAME( "output.csv" ){} //!< Automatic constructor
    std::string input_r1_fname; //!< Filename for forward reads, can be a .zip archive or a regular fastq file.
    std::string input_r2_fname; //!< Filename for reverse reads.
    std::string library_fname; //!< Filename containing a FASTA file containing a library of amino acid peptide sequences.
    std::string output_fname; //!< Filename of file to write output to.
    std::string samplelist_fname; //!< Name of tab-delimited file containing a list of samples.
    std::size_t seq_start; //!< Start index of where the sequences are expected to be found within the reads.
    std::size_t seq_len; //!< The length (how many nucleotides) of the designed peptides. 
    std::string index_fname; //!< Name of file containing indexed sequences.
    std::size_t max_mismatches; //!< The maximum number of mismatches to tolerate when mapping imperfectly-read reads.
    /**
     * The number of fastq records to read per loop. A higher value here will result in higher memory usage by the program.
     * However, higher values can also result in better performance as fewer disk accesses are performed.
     **/
    long int read_per_loop; 
    const long int DEFAULT_READ_PER_LOOP;
    const std::string DEFAULT_OUTPUT_FNAME;

    /**
     * Returns the arguments that are stored by 
     * the options object.
     **/
    std::string get_arguments();
};

#endif //OPTIONS_HH_INCLUDED
