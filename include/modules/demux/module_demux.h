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
     * @param input_map unordered_map which will have sequences as keys, and vectors of size_t as values. 
     * @param seqs The sequences that will be added to the unordered_map
     * @param num_samples The number of samples that were processed in the sequencing run.

     **/
    void add_seqs_to_map( parallel_map<sequence, std::vector<std::size_t>>& input_map, std::vector<sequence>& seqs, size_t num_samples );

};
#endif // MODULE_DEMUX_HH_INCLUDED
