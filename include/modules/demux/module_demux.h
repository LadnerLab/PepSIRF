#ifndef MODULE_DEMUX_HH_INCLUDED
#define MODULE_DEMUX_HH_INCLUDED
#include <vector>
#include <string>
#include <fstream>
#include <algorithm>
#include <omp.h>
#include <parallel_hashmap/phmap.h>

#include "options_demux.h"
#include "module.h"
#include "sequence.h"
#include "fasta_parser.h"
#include "fastq_parser.h"
#include "time_keep.h"

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

};
#endif // MODULE_DEMUX_HH_INCLUDED
