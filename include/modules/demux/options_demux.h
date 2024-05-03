#ifndef OPTIONS_DEMUX_HH_INCLUDED
#define OPTIONS_DEMUX_HH_INCLUDED
#include <string>
#include <tuple>
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
    std::string library_fname; //!< Filename containing a FASTA file containing a librry of amino acid peptide sequences.
    std::string output_fname; //!< Filename of file to write output to.
    std::string aggregate_fname; //!< Filename to write aggregate counts to.
    std::string flexible_idx_fname;
    bool pos_toggle;
    std::string samplelist_fname; //!< Name of tab-delimited file containing a list of samples.
    std::string samplename;
    std::string indexes;
    std::vector<std::string> sample_indexes;
    std::string diagnostic_fname;
    std::tuple<std::size_t, std::size_t, std::size_t> index1_data; //!< 0 = start, 1 = len, 2 = num_mismatches
    std::tuple<std::size_t, std::size_t, std::size_t> index2_data; //!< 0 = start, 1 = len, 2 = num_mismatches
    std::tuple<std::size_t, std::size_t, std::size_t> seq_data;     //!< 0 = start, 1 = len, 2 = num_mismatches
    std::string index_fname; //!< Name of file containing indexed sequences.
    std::string concatemer; //!< Concatenated primer sequences, we look for this in our reads to determine whether a peptide exists
    int phred_base;
    int min_phred_score;
    int num_indexes;
    std::string replicate_info_fname;

    bool translation_aggregation;
    std::string fastq_out;

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

    /**
     * Set the info of one of the f_index, r_index or seq data members.
     * Parses three comma-separated values from the info string, stores them in
     * the designated member.
     * @param member Pointer to either options_demux::f_index_data,
     *        options_demux::r_index_data, or options_demux::seq_data.
     *        The first item in the tuple pointed to by member will
     *        store the start-index of the specified item, the second the
     *        length and the third will store the number of mismatches
     *        to tolerate for this item.
     * @param info The string to parse. Must be 3 comma-separated
     *        integer values.
     **/
    void set_info( std::tuple<std::size_t, std::size_t, std::size_t>
                   options_demux:: *member, std::string info
                 );

    std::string tup_to_string( std::tuple<std::size_t,
                               std::size_t,
                               std::size_t>& data
                             );

};

#endif //OPTIONS_HH_INCLUDED
