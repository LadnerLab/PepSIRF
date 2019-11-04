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

    std::size_t id_index;

    /**
     * Name of the file to write output to.
     **/
    std::string output_fname;

    /**
     * Flag to penalize kmers based upon their frequency in the 
     * designed peptides. 
     **/
    bool penalize_kmers;

    /**
     * Name of file that contains names of enriched peptides,
     * one per line.
     **/
    std::string enriched_fname;

    /**
     * Threshold value for a peptide to not be filtered out.
     **/
    std::size_t threshold;

    /**
     * Directory name to write round scores/counts to
     **/
    std::string orig_scores_dname;

    /**
     * Score threshold to determine 
     * when two species are tied. 
     **/
    double score_tie_threshold;

    /**
     * Threshold to determine when 
     * two species have 'significant' 
     * overlap. Note that significant is not 
     * in the statistical sense, but is defined 
     * by either number or percentage of shared 
     * peptides.
     **/
    double score_overlap_threshold;

    /**
     * Name of the file to write a 
     * map of what peptides were assigned 
     * to what species to.
     **/
    std::string species_peptides_out;

    /**
     * If this value is true only one thread will be used for 
     * operations. Otherwise, two will be used.
     **/
    bool single_threaded;

    /**
     * Flag whether to use fractional scoring.
     **/
    bool fractional_scoring;

    /**
     * Flag whether to use 
     * summation scoring.
     **/
    bool summation_scoring;

    /**
     * Flag designating whether to use
     * score filtering instead of 
     * count filtering.
     **/
    bool score_filtering;

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

};

#endif // OPTIONS_DECONV_HH_INCLUDED
