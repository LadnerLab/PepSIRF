#ifndef MODULE_DEMUX_HH_INCLUDED
#define MODULE_DEMUX_HH_INCLUDED
#include <vector>
#include <string>
#include <fstream>
#include <limits>
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
#include "sequence_indexer.h"
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
     * Output is a tab-separated file, one line per sequence 
     * where each tab-separated entry is the count for a certain 
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


    /**
     * Method to zero a vector of size_t elements.
     * @param vec Pointer to the vector to zero.
     * @pre vec must have been initialized with some number of 
     *      elements, each element from vec[ 0 ] to vec[ size - 1 ]
     *      will be zero'd out.
     **/
    void _zero_vector( std::vector<std::size_t>* vec );

    std::size_t _get_read_length( std::ifstream& ofile );
    std::string _create_origin( std::size_t read_length );

    /**
     * Find the sequence mapped to, allowing for either one shift to the left or right,
     * or up to and including num_mism substitutions. 
     *W e first check to see if the match is exact, no bases were 
     * added, deleted, or changed during synthesization and reading.
     * Then, we shift one to the left and check again, and one to the right and check again.
     * Finally, we check to see if a match is found at the original location, but up to and 
     * include num_mism substitutions occurred. 
     * If no match is found, we return map.end()
     * @param map A map containing sequences as the keys, and anything as the value. 
     * @param probe_seq The sequence we are looking for in the map.
     * @param idx A sequence_indexer to query if no exact match is found when 
     *        searching from either shift. Note that the fuzzy match is only searched for 
     *        in the original expected position of the sequence.
     * @param num_mism The maximum number of mismatches (i.e. the maximum hamming distance) 
     *        to tolerate when searching for imperfectly-matched sequences. Any matches that are 
     *        found within this distance are included, and the minimum distance is used.  
     * @param f_start The start index at which we should look for probe_seq.
     * @param f_len The expected length of probe_seq, we search map for the substring of probe_seq
     *        starting at position f_start and ending at position f_len.
     * @returns Iterator to the match if found, map.end() otherwise
     **/
    template<class M>
        typename M::iterator _find_with_shifted_mismatch( M &map,
                                                          sequence probe_seq,
                                                          sequence_indexer& idx, std::size_t num_mism,
                                                          std::size_t f_start, std::size_t f_len
                                                         )
        {
            std::vector<std::pair<sequence *, int>> query_matches;
            sequence *best_match = nullptr;
            sequence seq_temp;

            unsigned int num_matches = 0;

            std::string substr = probe_seq.seq.substr( f_start, f_len );

            // Note: hash( sequence& seq ) = hash( seq.seq )
            auto temp = map.find( sequence( "", substr ) );

            // first check for an exact match in the expected location
            if( temp != map.end() )
                {
                    return temp;
                }

            if( f_start > 0 ) // check that we are not shifting left from the beginning
                {
                    substr = probe_seq.seq.substr( f_start - 1, f_len );
                    temp = map.find( sequence( "", substr ) );

                    if( temp != map.end() )
                        {
                            return temp;
                        }
                }

            // shift one to the right, look for exact match.
            substr = probe_seq.seq.substr( f_start + 1, f_len );
            temp = map.find( sequence( "", substr ) );

            if( temp == map.end() )
                {
                    substr = probe_seq.seq.substr( f_start, f_len );
                    seq_temp = sequence( "", substr );
                    num_matches = idx.query( query_matches,
                                             seq_temp,
                                             num_mism
                                           );
                    if( num_matches )
                        {
                            best_match = _get_min_dist( query_matches );
                            temp = map.find( *best_match );
                        }
                }

            // return the match, if found. Otherwise, map.end() is returned
            return temp;
        }

    /**
     * Gets the minimum distance from a vectors, returns a pointer to the 
     * sequence who has the minimum distance.
     * @note Because a max heap is used, returns the first item in the query vector.
     * @returns The item with the smallest value in the vector.
     **/
    sequence *_get_min_dist( std::vector<std::pair<sequence *, int>>& matches );


};


namespace demux
{
    void create_index_map( sequential_map<sequence, sample>& map,
                           std::vector<sequence>& index_seqs,
                           std::vector<sample>& samplelist
                         );

}

#endif // MODULE_DEMUX_HH_INCLUDED
