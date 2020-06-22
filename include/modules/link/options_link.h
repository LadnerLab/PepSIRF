#ifndef OPTIONS_LINK_HH_INCLUDED
#define OPTIONS_LINK_HH_INCLUDED
#include "options.h"
#include <string>

class options_link : public options
{
 public:
    options_link();

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
     * The kmer size to use
     * when mapping peptides.
     **/
    std::size_t k;

    /**
     * The name of the file to
     * write the id name map to.
     **/
    std::string id_name_map_fname;

    std::string output_fname;

    std::string metadata_fname;

    std::size_t id_index;

    /**
     * Flag to penalize kmers based upon their frequency in the
     * designed peptides.
     **/
    bool penalize_kmers;

    std::size_t kmer_size;

};

#endif // OPTIONS_LINK_HH_INCLUDED
