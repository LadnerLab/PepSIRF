#ifndef MODULE_DEMUX_HH_INCLUDED
#define MODULE_DEMUX_HH_INCLUDED
#include <vector>
#include <string>
#include <fstream>
#include <limits>
#include <algorithm>
#include <unordered_map>
#include <boost/algorithm/string.hpp>
#include "omp_opt.h"

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
#include <iomanip>
#include "fif_parser.h"

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
     * Writes output to diagnostic_fname. Output as tab-delimited file, with three columns, samplename, index pair matches, variable region matches.
     * Output is optional, defaulted as unused.
    **/
    void write_diagnostic_output( options_demux* d_opts,
                                  phmap::parallel_flat_hash_map<sample, std::vector<std::size_t>>& diagnostic_map );

    void create_diagnostic_map( bool reference_dependent,
                                phmap::parallel_flat_hash_map<sample,std::vector<std::size_t>>& diagnostic_map,
                                std::vector<sample> samplelist );

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
     * @param duplicate_map map of dna tags. contains the number of each dna tag that 
     *        appears in a run. used to determine the samples included in the output.
     **/
    void write_outputs( std::string outfile_name,
                        parallel_map<sequence, std::vector<std::size_t>*>& seq_scores,
                        std::map<std::string, std::size_t> duplicate_map,
                        bool ref_dependent,
                        std::vector<sample>& samples
                      );

    /**
     * Writes outputs to the outfile_name
     * 
     */


    /**
     * Method to zero a vector of size_t elements.
     * @param vec Pointer to the vector to zero.
     * @pre vec must have been initialized with some number of
     *      elements, each element from vec[ 0 ] to vec[ size - 1 ]
     *      will be zero'd out.
     **/
    void _zero_vector( std::vector<std::size_t>* vec );

    /**
     * Find the sequence mapped to, allowing for either one shift to the right or up to two shifts to the left,
     * or up to and including num_mism substitutions.
     * We first check to see if the match is exact, no bases were
     * added, deleted, or changed during synthesization and reading.
     * Then, we shift one and then two to the left and check again, and one to the right and check again.
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
    typename M::iterator _find_with_shifted_mismatch( M& map,
                                                      sequence probe_seq,
                                                      sequence_indexer& idx, std::size_t num_mism,
                                                      std::size_t f_start, std::size_t f_len
                                                    )
        {
            std::vector<std::pair<sequence *, int>> query_matches;
            sequence *best_match = nullptr;
            sequence seq_temp;
            std::string empty_string = "";
            unsigned int num_matches = 0;

            std::string substr = probe_seq.seq.substr( f_start, f_len );

            // Note: hash( sequence& seq ) = hash( seq.seq )
            auto temp = map.find( sequence( empty_string, substr ) );

            // first check for an exact match in the expected location
            if( temp != map.end() )
                {
                    return temp;
                }

            if( f_start > 0 ) // check that we are not shifting left from the beginning
                {
                    substr = probe_seq.seq.substr( f_start - 1, f_len );
                    temp = map.find( sequence( empty_string, substr ) );

                    if( temp == map.end() )
                        {
                            substr = probe_seq.seq.substr( f_start - 2, f_len );
                            temp = map.find( sequence( empty_string, substr ) );
                        }

                    if( temp != map.end() )
                        {
                            return temp;
                        }
                }

            // shift one to the right, look for exact match.
            // but only if we have enough substring to search
            if( f_start + 1 + f_len <= probe_seq.seq.length() )
                {
                    substr = probe_seq.seq.substr( f_start + 1, f_len );
                    temp = map.find( sequence( empty_string, substr ) );
                }

            // look for a match at the expected coordinates within
            // the number of mismatches that are tolerated
            if( num_mism > 0 && temp == map.end() )
                {
                    substr = probe_seq.seq.substr( f_start, f_len );
                    seq_temp.set_name( empty_string );
                    seq_temp.set_seq( substr );
                    num_matches = idx.query( query_matches,
                                             seq_temp,
                                             num_mism
                                           );
                    if( num_matches
                        && !_multiple_best_matches( query_matches ) )
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
     * @note Because matches is sorted, returns the first item in the vector
     * @pre  matches must be a sorted vector, sorted in non-decreasing order on the second item
     *       in each pair
     * @returns The item with the smallest value in the vector.
     **/
    sequence *_get_min_dist( std::vector<std::pair<sequence *, int>>& matches );

    /**
     * @note Creates two unordered multimaps based on DNA tag sequences referenced by the samplelist index columns.
     * The 'seq_lookup' map will contain the associated sequence to the ids in the barcodes/DNA tags file (--index-column).
     * Note only ids referenced in the samplelist file will be included.
     * 'map' contains concatenated DNA tag sequences, order specified by the index sample column name (--fif or --sIndex).
     * Example. The concatenation would follow as: Index1 sequence + Index2 sequence + Index3 sequence and so on. The accumulation of the sequences
     * is added to 'map' from the first index to the last index id for each sample. So three elements: sequence Index1, sequence Index1 + Index2,
     * and sequence Index1 + Index2 + Index3 will all be added.
     * @param map unordered multimap with unsorted sequence and samples. Samples organized into buckets where key sequence is identical.
     * This allows access to individual elements directly by the sequence object.
     * @param seq_lookup identical in concept to the map.
     **/
    void create_index_map( sequential_map<sequence, sample>& map,
                           std::vector<sequence>& dna_tags,
                           std::vector<sample>& samplelist,
                           sequential_map<sequence, sample>& seq_lookup
                         );

    /**
     * Determine if a sequence has had more than one best match. We say for sequence a
     * that a has multiple best matches iff the minimum of the second item in each
     * pair is not unique.
     * @param matches A vector of pairs where the first item is a pointer to sequence,
     *        and the second an integer count.
     * @returns boolean true if the minimum in matches appears more than once,
     *          false
     * @pre matches should be sorted
     **/
    bool _multiple_best_matches( std::vector<std::pair<sequence *, int>>& matches );

    /**
     * Aggregate output counts, creating count data at the aa-level.
     * When multiple encodings are created for each nt sequence, output
     * is reported at the aa-level. By summing the counts for each encoding of
     * each aa in the design, we get counts at the aa-level.
     * @param agg_map Reference to map where aggregate counts will be stored.
     * @param count_map Reference to map that contains counts at the nt-level.
     * @param num_samples The number of samples for which the design has been
     *        sequenced. This is equivalent to the number of lines in the samplelist
     *        file.
     * @pre The names of keys in count_map are formatted as such:
     *      NAME-ID, where NAME is the name of the aa sequence, and ID is the
     *      ID of the particular encoding of that sequence. Note the '-' separating
     *      the two identifiers.
     **/
    void aggregate_counts( parallel_map<sequence, std::vector<std::size_t>*>& agg_map,
                           parallel_map<sequence, std::vector<std::size_t>*>& count_map,
                           std::size_t num_samples
                         );

    /**
     * Write FASTQ output, appends fastq data to file
     * 
     * @param sequence Reference to fastq sequence for writing
     * @param sample_name sample file name to append
     */

    void write_fastq_output( std::map<std::string, std::vector<fastq_sequence>> samp_map,
                             std::string outfile
                           )
        {
            for( auto samp : samp_map ) 
                {
                    std::string outfile_path = "";
                    outfile_path = outfile + "/" + samp.first + ".fastq";
                    std::ofstream output;
                    output.open(outfile_path, std::ios_base::app);

                    for( auto fastq_seq : samp.second )
                        {
                            output << fastq_seq.name << "\n";
                            output << fastq_seq.seq << "\n";
                            output << "+" << "\n";
                            output << fastq_seq.scores << "\n";
                        }

                    output.close();
                }
        }


    /**
     * Aggregate output counts, creating count data at the aa-level.
     * When multiple encodings are created for each nt sequence, output
     * is reported at the aa-level. By summing the counts for each encoding of
     * each aa in the design, we get counts at the aa-level.
     * @param agg_map Reference to map where aggregate counts will be stored.
     * @param count_map Reference to map that contains counts at the nt-level.
     * @param num_samples The number of samples for which the design has been
     *        sequenced. This is equivalent to the number of lines in the samplelist
     *        file.
     * @note This method renames sequences to have the name of the sequences that
     *       resulted from the translation.
     **/
    template<typename Translator,
             typename Map
        >
        void aggregate_translated_counts( const Translator& translator,
                                          Map& agg_map,
                                          Map& count_map,
                                          size_t num_samples
                                        )
        {
            sequence current;
            std::vector<std::size_t> *current_vec_agg   = nullptr;
            std::vector<std::size_t> *current_vec_count = nullptr;
            sequential_map<std::string, std::vector<std::size_t>*> ptr_map;

            for( auto iter = count_map.begin(); iter != count_map.end(); ++iter )
                {

                    try
                        {
                            #pragma GCC diagnostic push
                            #pragma GCC diagnostic ignored "-Wdeprecated-copy"
                            current = translator( iter->first );
                            current.name = current.seq;
                            #pragma GCC diagnostic pop
                        }
                    catch( std::out_of_range& e )
                        {
                            continue;
                        }

                    // if the new name not in agg_map
                    if( ptr_map.find( current.name ) == ptr_map.end() )
                        {
                            // add the sequence to agg_map, initialize its data
                            agg_map[ current ] = new std::vector<std::size_t>( num_samples );
                            ptr_map[ current.name ] = agg_map[ current ];
                            std::fill( agg_map[ current ]->begin(), agg_map[ current ]->end(), 0 );
                        }

                    current_vec_agg   = ptr_map[ current.name ];
                    current_vec_count = iter->second;

                    // add the counts from the untrimmed count_map entry to
                    // the trimmed agg_map entry
                    for( std::size_t index = 0; index < num_samples; ++index )
                        {
                            current_vec_agg->at( index ) += current_vec_count->at( index );
                        }
                }


        }

};


#endif // MODULE_DEMUX_HH_INCLUDED
