#ifndef OPTIONS_LINK2_HH_INCLUDED
#define OPTIONS_LINK2_HH_INCLUDED
#include "options.h"
#include <string>

class options_link2 : public options
{
 public:
    options_link2();

    std::string get_arguments();

    /**
     * Name of the fasta file containing
     * protein sequences
     **/
    std::string prot_file_fname;

    /**
     * Name of the fasta file containing peptide
     * sequences.
     **/
    std::string peptide_file_fname;

    /**
     * Metadata info: 1) name of tab-delimited metadata file, 
     * 2) header for column containing protein sequence name and 
     * 3) header for column containing ID to be used in creating the linkage map
     */
     std::string metadata_info;

    /**
     * The kmer size to use
     * when mapping peptides.
     **/
    std::size_t kmer_size;

    /**
     * Number of positions across which the match can span.
     **/
    std::size_t span;

    std::string output_fname;
};

#endif // OPTIONS_LINK2_HH_INCLUDED
