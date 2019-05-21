#ifndef OPTIONS_DEMUX_HH_INCLUDED
#define OPTIONS_DEMUX_HH_INCLUDED
#include <string>
#include "options.h"

/*! Data class to contain and handle 
 * arguments that will be passed in from the 
 * command-line.
*/
class options_demux: public options
{
public:
	options_demux(); //!< Automatic constructor
    std::string input_r1_fname; //!< Filename for forward reads, can be a .zip archive or a regular fastq file.
    std::string input_r2_fname; //!< Filename for reverse reads, can be a .zip archive or regular fastq file.
    std::string library_fname; //!< Filename containing a FASTA file containing a library of amino acid peptide sequences.
};

#endif //OPTIONS_HH_INCLUDED
