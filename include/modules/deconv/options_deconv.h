#ifndef OPTIONS_DECONV_HH_INCLUDED
#define OPTIONS_DECONV_HH_INCLUDED

#include "logger.h"
#include "options.h"

#include <string>
#include <tuple>
#include <vector>
#include <boost/algorithm/string.hpp>

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
     * Suffix to add to files written to
     * output directory in batch mode.
     **/
    std::string outfile_suffix;

    /**
     * Suffix to add to map files.
     **/
    std::string map_suffix;

    /**
     * Name of the file to write output to.
     **/
    std::string output_fname;

    /**
     * Flag to penalize kmers based upon their frequency in the
     * designed peptides.
     **/
    bool penalize_kmers;

    bool remove_file_types;

    /**
     * Name of file that contains names of enriched peptides,
     * one per line.
     **/
    std::string enriched_fname;

    /**
     * Name of file that contains threshold values for each 
     * peptide species to not be filtered out.
     **/
    std::string thresholds_fname;

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
     * Flag stores scoring strategy. Default: summation.
     **/
    std::string scoring_strategy;

    /**
     * Flag designating whether to use
     * score filtering instead of
     * count filtering.
     **/
    bool score_filtering;

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
     * Tuple containing information for parsing user-defined ID name map
     */
    std::tuple<std::string, std::string, std::string>
        custom_id_name_map_info; //!< 0 = ID name map, 1 = Taxon ID column name, 2 = Sequence column name

    /**
     * Name of file containing NCBI taxID to taxon name mappings
     */
    std::string id_name_map_fname;

    /**
     * The expected ending of the file name
     * with enriched peptides
     **/
    std::string enriched_file_ending;

    /**
     * Translates information from string 'info' into a tuple of length 3 and
     * any type.
     * Note: is probably only going to be used to move custom ID name map
     * information.
     * @param tup Reference to tuple which will be filled
     * @param info Comma-delimited string of information
     */
    // TODO: make into a generic function in options class,
    // or a virtual method - whichever is more practical for the future
    void set_info(
        std::tuple<std::string, std::string, std::string>& tup,
        std::string info
    );

    /**
     */
    // TODO: make into a generic function in options class,
    // or a virtual method - whichever is more practical for the future
    std::string tuple_to_string(std::tuple<std::string, std::string, std::string>& tup);
};

#endif /* OPTIONS_DECONV_HH_INCLUDED */

