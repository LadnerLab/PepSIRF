#ifndef MODULE_DEMUX_HH_INCLUDED
#define MODULE_DEMUX_HH_INCLUDED
#include <vector>
#include <string>
#include <fstream>
#include <algorithm>
#include <omp.h>

#include "options_demux.h"
#include "maps.h"
#include "module.h"
#include "sequence.h"
#include "fasta_parser.h"
#include "fastq_parser.h"
#include "time_keep.h"
#include "samplelist_parser.h"
#include "sample.h"


/**
 * Class for running the demultiplex module. Given a file of reads and a file containing 
 * a designed set of oligos, maps reads to a reference and counts how many times each reference 
 * sequence was read. 
 **/
class module_demux : public module
{
 public:

    std::string name;

    /**
     * Default constructor, sets the 'name' member of the 
     * class to 'demux'.
     **/
    module_demux();

    /**
     * Runs the demux (demultiplex) module. 
     * @param opts Pointer to an instance of the 
     *        'options_demux' class whose values have been 
     *        initialized. 
     * @pre Opts is an instance of the 'options_demux' class 
     *      that has been initialized. If opts is unitialized, 
     *      undefined behavior will result. 
     **/
    void run( options *opts );


    std::string get_name();

    /**
     * Adds sequences to an unordered map, where the key is the string sequence, and the value is a vector of 
     * size_t counts. Each count will be used to keep track of how many times each sequence is found per sample.
     * @param input_map parallel_map which will have sequences as keys, and vectors of size_t as values. 
     * @param seqs The sequences that will be added to the unordered_map
     * @param num_samples The number of samples that were processed in the sequencing run.

     **/
    void add_seqs_to_map( parallel_map<sequence, std::vector<std::size_t>*>& input_map, std::vector<sequence>& seqs, size_t num_samples );


    /**
     * Writes output to the outfile_name.
     * Output is a comma-separated file, one line per sequence 
     * where each comma-separated entry is the count for a certain 
     * sample. 
     * @note A header is written to the file labelling each column.
     * @param outfile_name Name of file to write to.
     * @param seq_scores A map that couples sequences with a vector of their scores, where 
     *        the i'th entry in the seq_scores for sequence j is the number of times that 
     *        sequence j was found in sample i.
     * @param samples vector of samples. This vector is used to label each of the samples in the 
     *        vector of each seq_score. Note that samples[ i ].id must equal j[ i ] for each 
     *        j = 1, 2, ... j.size(), i.e. The id of a sample must correspond with its entry in 
     *        the count vector.
     **/
    void write_outputs( std::string outfile_name,
                        parallel_map<sequence, std::vector<std::size_t>*>& seq_scores,
                        std::vector<sample>& samples
                      );


};

namespace demux
{
    void create_index_map( sequential_map<std::string, sample>& map,
                           std::vector<sequence>& index_seqs,
                           std::vector<sample>& samplelist
                         );

}

#endif // MODULE_DEMUX_HH_INCLUDED
