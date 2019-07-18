#ifndef OPTIONS_DECONV_HH_INCLUDED
#define OPTIONS_DECONV_HH_INCLUDED

#include "options.h"

#include <string>
#include <vector>

class options_deconv : public options
{
 public:
    
    /**
     * Returns the arguments that are stored by 
     * the options object.
     **/
    std::string get_arguments();

    /**
     * Name of the file containing lines linking 
     * peptides to the species each peptide shares 
     * a 7-mer with.
     **/
    std::string linked_fname;

    /**
     * Name of the file to write output to.
     **/
    std::string output_fname;

    /**
     * Name of file that contains names of enriched peptides,
     * one per line.
     **/
    std::string enriched_fname;

    /**
     * Threshold value for a peptide to not be filtered out.
     **/
    float threshold;

    /**
     * If this value is true only one thread will be used for 
     * operations. Otherwise, two will be used.
     **/
    bool single_threaded;

    /**
     * Flag whether to use fractional scoring.
     **/
    bool fractional_scoring;

    bool summation_scoring;

    /**
     * Bool to determine whether to 
     * create linkage.
     **/
    bool create_linkage;

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

    std::size_t k;

    std::string id_name_map_fname;

};

#endif // OPTIONS_DECONV_HH_INCLUDED
