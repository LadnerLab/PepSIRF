#include <sstream>
#include <stdexcept>

#include "file_io.h"
#include "module_demux.h"
#include "fastq_sequence.h"
#include "fastq_score.h"
#include "et_search.h"
#include "nt_aa_translator.h"

module_demux::module_demux()
{
    name = "Demux";
}

void module_demux::run( options *opts )
{
    // options from the command line
    options_demux *d_opts = (options_demux*) opts;

    std::string index_str;

    const std::size_t forward_start  = std::get<0>( d_opts->f_index_data );
    const std::size_t forward_length = std::get<1>( d_opts->f_index_data );

    const std::size_t reverse_start  = std::get<0>( d_opts->r_index_data );
    const std::size_t reverse_length = std::get<1>( d_opts->r_index_data );


    std::size_t read_index = 0;

    struct time_keep::timer total_time;
    parallel_map<sequence, std::vector<std::size_t>*> reference_counts;
    parallel_map<sequence, std::size_t> non_perfect_match_seqs;

    sequential_map<sequence, sample> index_map;

    omp_set_num_threads( opts->num_threads );

    total_time.start();;

    // vector to store the .fna sequences that represent a designed library
    std::vector<sequence> library_seqs;
    std::vector<sequence> index_seqs;
    std::vector<fastq_sequence> r2_seqs;
    std::vector<fastq_sequence> reads;

    reads.reserve( d_opts->read_per_loop );

    // create parsers for the encoded library and the fastq reads.
    fasta_parser fasta_p;
    fastq_parser fastq_p;

    samplelist_parser samplelist_p;
    std::vector<sample> samplelist = samplelist_p.parse( d_opts->samplelist_fname );

    // open the read file so we can iterate over it
    std::ifstream reads_file( d_opts->input_r1_fname, std::ios_base::in );
    std::ifstream r2_reads;

    bool r2_reads_included = d_opts->input_r2_fname.length() > 0;

    if( r2_reads_included )
        {
            r2_reads.open( d_opts->input_r2_fname, std::ios_base::in );
        }

    const bool reference_dependent = !d_opts->library_fname.empty();

    if( reference_dependent )
        {
            library_seqs   = fasta_p.parse( d_opts->library_fname );
        }
    else
        {
            std::cout << "WARNING: A set of reference sequences was not provided. "
                      << "Demux will be run in reference-independent mode, where each "
                      << "read is treated as its own reference.\n";
        }

    index_seqs     = fasta_p.parse( d_opts->index_fname );

    create_index_map( index_map, index_seqs, samplelist );
    sequential_map<sequence, sample> sample_table;
    sequential_map<sequence,sample>::size_type num_samples = index_map.size();

    add_seqs_to_map( reference_counts, library_seqs, samplelist.size() );

    std::string adapter;
    std::string nuc_seq;

    parallel_map<sequence, std::vector<std::size_t>*>::iterator seq_iter;

    std::size_t seq_start  = std::get<0>( d_opts->seq_data );
    std::size_t seq_length = std::get<1>( d_opts->seq_data );

    std::size_t processed_total   = 0;
    std::size_t processed_success = 0;
    std::size_t concatemer_found  = 0;

    sequence_indexer lib_idx;
    sequence_indexer index_idx;

    // index the forward index sequences and the library
    // sequences
    lib_idx.index( library_seqs );
    index_idx.index( index_seqs );

    std::size_t sample_id = 0;
    sequential_map<sequence, sample>::iterator r_idx_match;
    sequential_map<sequence, sample>::iterator f_idx_match;

    std::istream *reads_ptr = &reads_file;
    std::istream *r2_reads_ptr = &r2_reads;

    #ifdef ZLIB_ENABLED

    pepsirf_io::gzip_reader gzip_r1_reader( reads_file );
    pepsirf_io::gzip_reader gzip_r2_reader( r2_reads );


    #endif

    if( pepsirf_io::is_gzipped( reads_file ) )
        {

    #ifdef ZLIB_ENABLED
            reads_ptr = &gzip_r1_reader;
    #else
            throw std::runtime_error( "A gzipped file was supplied, but "
                                      "boost zlib is not available. "
                                      "Please either configure boost with zlib "
                                      "or unzip the file before attempting to demux."
                                    );
    #endif
        }

    if( r2_reads_included
        && pepsirf_io::is_gzipped( r2_reads )
      )
        {

    #ifdef ZLIB_ENABLED
            r2_reads_ptr = &gzip_r2_reader;
    #else
            throw std::runtime_error( "A gzipped file was required, but "
                                      "boost zlib is not available. "
                                      "Please either configure boost with zlib "
                                      "or unzip the file before attempting to demux."
                                    );
    #endif
        }


    std::istream& reads_file_ref = *reads_ptr;
    std::istream& r2_reads_ref   = *r2_reads_ptr;
    std::size_t f_idx_match_count = 0;
    std::size_t both_idx_match_count = 0;
    std::size_t var_reg_match_count = 0;
    d_opts->total_pair_matches = {0};
    d_opts->var_region_matches = {0};
    // while we still have reads to process, loop terminates when no fastq records
    // are read from the file
    while( fastq_p.parse( reads_file_ref, reads, d_opts->read_per_loop  ) )
        {

            if( r2_reads_included )
                {
                    // get the index sequences for this set of sequences
                    fastq_p.parse( r2_reads_ref, r2_seqs, d_opts->read_per_loop );
                }

           #pragma omp parallel for private( seq_iter, nuc_seq, read_index, index_str, adapter, sample_id,  \
                                              r_idx_match, f_idx_match ) \
               shared( seq_start, seq_length, d_opts, reference_counts, library_seqs, index_seqs, r2_seqs ) \
                reduction( +:processed_total, processed_success, concatemer_found ) schedule( dynamic )
            for( read_index = 0; read_index < reads.size(); ++read_index )
                {
                    // lambda to return whether a match was found
                    // if no reverse indexes are used, then we only
                    // check the forward index match, otherwise
                    // we need to check both
                    auto match_found = [&]() -> bool
                        {
                            return
                            ( reverse_length == 0
                              && f_idx_match != index_map.end()
                              )
                            ||
                            ( reverse_length != 0
                                 && f_idx_match != index_map.end()
                                 && r_idx_match != index_map.end()
                               );


                        };


                    // checks to see that a concatemer was included
                    auto found_concatemer = [&]() -> bool
                        {
                            return
                                     d_opts->concatemer.length() > 0
                                       && reads[ read_index ].seq.find( d_opts->concatemer,
                                                                        0
                                                                        ) != std::string::npos;

                        };

                    auto quality_match = [&]() -> bool
                        {
                            return ( !d_opts->min_phred_score
                                     || ( fastq_score::get_avg_score( reads[ read_index ]
                                                                      .scores.begin() + seq_start,

                                                                      reads[ read_index ]
                                                                      .scores.begin() + seq_start + seq_length,
                                                                      d_opts->phred_base
                                                                    ) >= d_opts->min_phred_score
                                        )
                                   );
                        };

                    if( reverse_length > 0 )
                        {
                            if( r2_seqs.size() == 0 )
                                {
                                    r_idx_match = _find_with_shifted_mismatch( index_map, reads[ read_index ],
                                                                               index_idx, std::get<2>( d_opts->r_index_data ),
                                                                               reverse_start, reverse_length
                                                                             );
                            }
                            else
                                {
                                    r_idx_match = _find_with_shifted_mismatch( index_map, r2_seqs[ read_index ],
                                                                               index_idx, std::get<2>( d_opts->r_index_data ),
                                                                               reverse_start, reverse_length
                                                                             );
                                }
                        }

                    // get the forward and the reverse indexes from the sequence, grab the id
                    // for this index. This gives the index of the location at the sequence to increment
                    f_idx_match = _find_with_shifted_mismatch( index_map, reads[ read_index ],
                                                               index_idx, std::get<2>( d_opts->f_index_data ),
                                                               forward_start, forward_length
                                                             );
                    if( match_found()
                        && quality_match()
                      )
                        {
                            if( reverse_length == 0 )
                                {
                                    f_idx_match_count++;
                                }
                            else
                                {
                                    both_idx_match_count++;
                                }

                            using seq_map = parallel_map<sequence, std::vector<std::size_t>*>;

                            if( reference_dependent )
                                {
                                    et_seq_search<seq_map,true> library_searcher( lib_idx, reference_counts, num_samples );


                                    auto seq_match = library_searcher.find( reads[ read_index ],
                                                                            std::get<2>( d_opts->seq_data ),
                                                                            seq_start,
                                                                            seq_length
                                                                          );
                                    if( reverse_length == 0
                                        && seq_match != reference_counts.end()
                                        )
                                        {
                                            sample_id = f_idx_match->second.id;
                                            seq_match->second->at( sample_id )++;
                                            ++processed_success;
                                            var_reg_match_count++;
                                        }
                                    else if( reverse_length != 0
                                             && seq_match != reference_counts.end()
                                             )
                                        {
                                            std::string concat_idx = f_idx_match->first.seq
                                                + r_idx_match->first.seq;

                                            auto d_id = index_map.find( sequence( "", concat_idx ) );

                                            if( d_id != index_map.end() )
                                                {
                                                    sample_id = d_id->second.id;
                                                    seq_match->second->at( sample_id )++;
                                                    ++processed_success;
                                                    var_reg_match_count++;
                                                }
                                        }
                                    else if( seq_match == reference_counts.end()
                                             && found_concatemer() )
                                        {
                                            ++concatemer_found;
                                        }
                                }
                            else
                                {
                                    // increase count by one for every index match found
                                    d_opts->total_pair_matches[ read_index ]++;

                                    et_seq_search<seq_map,false> library_searcher( lib_idx, reference_counts, num_samples );


                                    if( reverse_length == 0 )
                                        {
                                            auto seq_match = library_searcher.find( reads[ read_index ],
                                                                                    std::get<2>( d_opts->seq_data ),
                                                                                    seq_start,
                                                                                    seq_length
                                                                                    );

                                            if( seq_match != reference_counts.end() )
                                                {
                                                    sample_id = f_idx_match->second.id;
                                                    seq_match->second->at( sample_id )++;
                                                    ++processed_success;
                                                    var_reg_match_count++;
                                                }
                                        }
                                    else if( reverse_length != 0 )
                                        {
                                            std::string concat_idx = f_idx_match->first.seq
                                                + r_idx_match->first.seq;

                                            auto d_id = index_map.find( sequence( "", concat_idx ) );

                                            if( d_id != index_map.end() )
                                                {
                                                    std::size_t n_found = reads[ read_index ].seq.find( "N" );
                                                    if( n_found == std::string::npos )
                                                        {
                                                            auto seq_match = library_searcher.find( reads[ read_index ],
                                                                                                    std::get<2>( d_opts->seq_data ),
                                                                                                    seq_start,
                                                                                                    seq_length
                                                                                                    );
                                                            if( seq_match != reference_counts.end() )
                                                                {

                                                                    sample_id = d_id->second.id;
                                                                    seq_match->second->at( sample_id )++;
                                                                    ++processed_success;
                                                                    var_reg_match_count++;
                                                                }
                                                        }
                                                }
                                        }
                                    else if( found_concatemer() )
                                        {
                                            ++concatemer_found;
                                        }

                                }

                        }
                    // record the number of records that are processed
                    ++processed_total;
                }

            reads.clear();
            r2_seqs.clear();

        }

    total_time.stop();

    std::cout << "Processed " << processed_total << " records in " << time_keep::get_elapsed( total_time ) << " seconds.\n";
    std::cout << processed_success << " records were found to be a match out of "
              << processed_total << " (" << ( (long double) processed_success / (long double) processed_total ) * 100
              << "%) successful.\n";
    std::cout << ( ( long double ) f_idx_match_count / ( long double ) processed_total ) * 100.00 << "% of reads from index 1 were found to be a match out of total.\n"
              << ( ( long double ) both_idx_match_count / ( long double ) processed_total ) * 100.00 << "% of reads from index 1 and index 2 were found to be a match out of total.\n"
              << ( ( long double ) var_reg_match_count / ( long double ) processed_total ) * 100.00 << "% of reads from variable sequence regions were found to be a match out of total.\n";

    if( d_opts->concatemer.length() > 0 )
        {
            std::cout << "The concatemer sequence was found " << concatemer_found << " times (" <<
                ( (long double) concatemer_found / (long double) processed_total ) * 100 << "% of total).\n";
        }

    if( d_opts->aggregate_fname.length() > 0  )
        {
            parallel_map<sequence, std::vector<std::size_t>*> agg_map;

            if( d_opts->translation_aggregation )
                {
                    nt_aa_translator<codon_aa_mappings::default_map_type>
                        translator( codon_aa_mappings::default_codon_aa_map );

                    aggregate_translated_counts( translator,
                                                 agg_map,
                                                 reference_counts,
                                                 samplelist.size()
                                               );

                }
            else
                {
                    aggregate_counts( agg_map, reference_counts, samplelist.size() );
                }

            write_outputs( d_opts->aggregate_fname,
                           agg_map,
                           samplelist
                         );
        }
    write_outputs( d_opts->output_fname, reference_counts, samplelist );
}


std::string module_demux::get_name()
{
    return name;
}

void module_demux::aggregate_counts( parallel_map<sequence, std::vector<std::size_t>*>& agg_map,
                                     parallel_map<sequence, std::vector<std::size_t>*>& count_map,
                                     size_t num_samples
                                   )
{
    auto iter = count_map.begin();
    std::size_t index = 0;
    std::vector<std::string> strs;
    sequence current;
    std::vector<std::size_t> *current_vec_agg   = nullptr;
    std::vector<std::size_t> *current_vec_count = nullptr;
    sequential_map<std::string, std::vector<std::size_t>*> ptr_map;

    constexpr int NUM_DELIMITED_ITEMS = 2;

    // for each in count map
    for( ; iter != count_map.end(); ++iter )
        {
        // trim the '-ID' from the name of the sequence
            boost::split( strs, iter->first.name, boost::is_any_of( "-") );

            std::string cpy = strs[ 0 ];
            current = sequence( cpy, iter->first.seq );

            // we are not going to get correct results
            // if names are not in the expected
            // format.
            if( strs.size() != NUM_DELIMITED_ITEMS )
                {

                    std::ostringstream err_msg;
                    err_msg << "The following nucleotide sequence is not "
                            << "formatted correctly in order to retrieve "
                            << "aa-level counts: \n"
                            << iter->first.name
                            << "\n";
                   throw std::runtime_error( err_msg.str() );
                }

            // if the trimmed name not in agg_map
            if( ptr_map.find( strs[ 0 ] ) == ptr_map.end() )
                {
                    // add the sequence to agg_map, initialize its data
                    agg_map[ current ] = new std::vector<std::size_t>( num_samples );
                    ptr_map[ strs[ 0 ] ] = agg_map[ current ];
                    _zero_vector( agg_map[ current ] );
                }

            current_vec_agg   = ptr_map[ current.name ];
            current_vec_count = iter->second;

            // add the counts from the untrimmed count_map entry to
            // the trimmed agg_map entry
            for( index = 0; index < num_samples; ++index )
                {
                    current_vec_agg->at( index ) += current_vec_count->at( index );
                }
        }

}
void module_demux::add_seqs_to_map( parallel_map<sequence, std::vector<std::size_t>*>& input_map, std::vector<sequence>& seqs, size_t num_samples )
{
    std::size_t index        = 0;

    input_map.reserve( seqs.size() );

    // #pragma omp parallel for private( index ) shared ( seqs, input_map, num_samples )
    for( index = 0; index < seqs.size(); ++index )
        {
            input_map[ seqs[ index ] ] = new std::vector<std::size_t>( num_samples );
        }
    for( auto& x : input_map )
        {
            std::fill( x.second->begin(),
                       x.second->end(),
                       0
                     );
        }
}

void module_demux::write_diagnostic_output( std::string outfile_name, options_demux* d_opts, std::vector<sample>& samples )
{

    std::ofstream outfile( outfile_name, std::ofstream::out );
    std::string line = "Sequence name\tRead Matches\tVariable Region Matches\n";
    outfile << line;
    std::size_t id_index = 0;
    for( const auto& sample_name : samples )
        {
            outfile << sample_name.name << "\t" << d_opts->total_pair_matches[ id_index ] << "\t" << d_opts->var_region_matches[ id_index ] << "\n";
            id_index++;
        }
    outfile.close();

}

void module_demux::write_outputs( std::string outfile_name,
                                  parallel_map<sequence, std::vector<std::size_t>*>& seq_scores,
                                  std::vector<sample>& samples
                                )
{
    std::ofstream outfile( outfile_name, std::ofstream::out );
    parallel_map<sequence, std::vector<std::size_t>*>::iterator seq_iter = seq_scores.begin();

    const std::string DELIMITER = "\t";
    const std::string NEWLINE   = "\n";


    std::size_t index        = 0;
    std::size_t second_index = 0;

    std::string header       = "Sequence name";
    header.append( DELIMITER );

    for( index = 0; index < samples.size() - 1; ++index )
        {
            header.append( samples[ index ].name );
            header.append( DELIMITER );
        }
    header.append( samples[ index ].name );
    header.append( NEWLINE );

    outfile << header;

    for( index = 0; index < seq_scores.size(); ++index )
        {
            const sequence& curr = seq_iter->first;
            const std::vector<std::size_t> *curr_counts = seq_iter->second;


            outfile << curr.name << DELIMITER;

            for( second_index = 0; second_index < samples.size() - 1; ++second_index )
                {
                    outfile << curr_counts->at( second_index ) << DELIMITER;
                }
            outfile << curr_counts->at( samples.size() - 1 ) << NEWLINE;

            ++seq_iter;
            delete curr_counts;
        }
    outfile.close();
}

void module_demux::_zero_vector( std::vector<std::size_t>* vec )
{
    std::size_t index = 0;
    for( index = 0; index < vec->size(); ++index )
        {
            vec->at( index ) = 0;
        }
}

void module_demux::create_index_map( sequential_map<sequence, sample>& map,
                                     std::vector<sequence>& index_seqs,
                                     std::vector<sample>& samplelist
                                   )
{
    sequential_map<std::string, std::string> st_table;
    std::string concat_seq = "";
    std::string seq1       = "";
    std::string seq2       = "";

    unsigned int index = 0;

    for( index = 0; index < index_seqs.size(); ++index )
        {
            st_table[ index_seqs[ index ].name ] = index_seqs[ index ].seq;
            map[ sequence( index_seqs[ index ].name,
                           index_seqs[ index ].seq
                         )
               ] = samplelist[ 0 ];
        }

    for( index = 0; index < samplelist.size(); ++index )
        {
            seq1 = st_table[ samplelist[ index ].get_first_id() ];
            seq2 = st_table[ samplelist[ index ].get_second_id() ];
            concat_seq = seq1 + seq2;

            map[ sequence( "", concat_seq ) ] = samplelist[ index ];
        }
}


bool module_demux::_multiple_best_matches( std::vector<std::pair<sequence *, int>>& matches )
{
    return matches.size() >= 2 &&
        std::get<1>( matches[ 0 ] ) == std::get<1>( matches[ 1 ] );
}

sequence *module_demux::_get_min_dist( std::vector<std::pair<sequence *, int>>& matches )
{
    return std::get<0>( *matches.begin() );
}

