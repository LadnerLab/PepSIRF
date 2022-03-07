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

    // fif use case
    std::vector<flex_idx> flexible_idx_data;
    std::vector<sequential_map<sequence, sample>::iterator> idx_match_list;

    if( d_opts->flexible_idx_fname.empty() )
        {
            flexible_idx_data.emplace_back( d_opts->sample_indexes[0],            // index name
                                            "r1",                                 // read direction
                                            std::get<0>( d_opts->index1_data ),   // start
                                            std::get<1>( d_opts->index1_data ),   // length
                                            std::get<2>( d_opts->index1_data ) ); // num mismatch
            flexible_idx_data.emplace_back( d_opts->sample_indexes[1],
                                            "r2",
                                            std::get<0>( d_opts->index2_data ),
                                            std::get<1>( d_opts->index2_data ),
                                            std::get<2>( d_opts->index2_data ));
        }
    else
        {
            fif_parser flex_parser;
            flexible_idx_data = flex_parser.parse( d_opts->flexible_idx_fname );
        }


    std::size_t read_index = 0;
    struct time_keep::timer total_time;
    parallel_map<sequence, std::vector<std::size_t>*> reference_counts;
    parallel_map<sequence, std::size_t> non_perfect_match_seqs;

    sequential_map<sequence, sample> index_map;
    //std::map<std::string, std::size_t> seq_duplicates;
    std::unordered_map<std::string, std::vector<std::tuple<std::string, std::string, std::size_t>>> diagnostic_map;

    omp_set_num_threads( opts->num_threads );

    total_time.start();

    // vector to store the .fna sequences that represent a designed library
    std::vector<sequence> library_seqs;
    std::vector<sequence> index_seqs;
    std::vector<fastq_sequence> r2_seqs;
    std::vector<fastq_sequence> reads;

    reads.reserve( d_opts->read_per_loop );
    // create parsers for the encoded library and the fastq reads.
    fasta_parser fasta_p;
    fastq_parser fastq_p;
    // parse samplelist
    samplelist_parser samplelist_p;

    if( !d_opts->flexible_idx_fname.empty() )
        {
            for( std::size_t curr_row = 0; curr_row < flexible_idx_data.size(); ++curr_row )
                {
                    d_opts->sample_indexes.emplace_back( flexible_idx_data[curr_row].idx_name );
                }
        }
    std::vector<sample> samplelist = samplelist_p.parse( d_opts );

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
    sequential_map<sequence, sample> seq_lookup;
    create_index_map( index_map, index_seqs, samplelist, seq_lookup);
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

    lib_idx.index( library_seqs );
    index_idx.index( index_seqs );

    std::size_t sample_id = 0;
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
    std::vector<std::pair<std::string,size_t>> index_match_totals;
    for( const auto& curr_index : flexible_idx_data )
        {
            index_match_totals.emplace_back( std::make_pair( curr_index.idx_name, 0 ) );
        }

    if( !d_opts->diagnostic_fname.empty() )
        {
            for( const auto sample : samplelist )
                {
                    std::vector<std::tuple<std::string, std::string, std::size_t>> barcodes;
                    for( std::size_t curr_barcode = 0; curr_barcode < sample.string_ids.size(); ++curr_barcode )
                        {
                            barcodes.emplace_back( std::make_tuple( index_match_totals[curr_barcode].first, sample.string_ids[curr_barcode], 0 ) );
                        }
                    if( reference_dependent )
                        {
                            barcodes.emplace_back( std::make_tuple("DNA Tag Matches", "", 0) );
                        }
                    diagnostic_map.emplace( sample.name, barcodes );
                }
        }

    while( fastq_p.parse( reads_file_ref, reads, d_opts->read_per_loop  ) )
        {

            if( r2_reads_included )
                {
                    fastq_p.parse( r2_reads_ref, r2_seqs, d_opts->read_per_loop );
                }

           #pragma omp parallel for private( seq_iter, nuc_seq, read_index, index_str, adapter, sample_id,  \
                                              idx_match_list ) \
               shared( seq_start, seq_length, d_opts, num_samples, reference_counts, library_seqs, index_seqs, r2_seqs) \
                reduction( +:processed_total, processed_success, concatemer_found ) \
                schedule( dynamic )

            for( read_index = 0; read_index < reads.size(); ++read_index )
                {
                    
                    auto match_found = [&]() -> bool
                        {
                            std::string concat_idx_seq = "";
                            
                            std::vector<std::string> matched_ids;
                            for( std::size_t curr_index = 0; curr_index < idx_match_list.size(); curr_index++ )
                                {
                                    // determine which indexes found matches
                                    sequential_map<sequence, sample>::iterator index_match = idx_match_list[curr_index];
                                    if( index_match != index_map.end() )
                                        {
                                            matched_ids.emplace_back( index_match->second.string_ids[curr_index] );
                                            concat_idx_seq.append( index_match->first.seq );
                                            if( index_map.find( sequence( "", concat_idx_seq ) ) != index_map.end() )
                                                {
                                                    index_match_totals[curr_index].second++;
                                                    // future: update_diagnostic_data - or something
                                                    if( !d_opts->diagnostic_fname.empty() )
                                                        {
                                                            for( auto& sample : diagnostic_map )
                                                                {
                                                                    std::size_t curr_barcode = 0;
                                                                    while( curr_barcode < matched_ids.size()
                                                                           && std::get<1>( sample.second[curr_barcode] ) == matched_ids[curr_barcode] )
                                                                        {
                                                                            // if the current samples barcode id matches the current id from the found match
                                                                            if( curr_barcode == curr_index )
                                                                                {
                                                                                    std::get<2>(sample.second[curr_barcode]) += 1;
                                                                                }
                                                                            ++curr_barcode;
                                                                        }
                                                                }
                                                        }
                                                }
                                            else
                                                {
                                                    return false;
                                                }
                                        }
                                    else
                                        {
                                            return false;
                                        }
                                }
                            return true;
                        };

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
                    for( const auto& index : flexible_idx_data )
                        {
                            sequential_map<sequence, sample>::iterator index_match;
                            if( index.read_name == "r1" )
                                {
                                    index_match = _find_with_shifted_mismatch( seq_lookup, reads[ read_index ],
                                                                            index_idx, index.num_mismatch,
                                                                            index.idx_start, index.idx_len
                                                                            );
                                }
                            else
                                {
                                    if( r2_seqs.size() == 0 )
                                        {
                                            index_match = _find_with_shifted_mismatch( seq_lookup, reads[ read_index ],
                                                                                    index_idx, index.num_mismatch,
                                                                                    index.idx_start, index.idx_len
                                                                                    );
                                        }
                                    else
                                        {
                                            index_match = _find_with_shifted_mismatch( seq_lookup, r2_seqs[ read_index ],
                                                                                    index_idx, index.num_mismatch,
                                                                                    index.idx_start, index.idx_len
                                                                                    );
                                        }
                                }

                            idx_match_list.push_back( index_match );
                        }

                    if( match_found()
                        && quality_match()
                      )
                        {
                            using seq_map = parallel_map<sequence, std::vector<std::size_t>*>;
                            sequential_map<sequence,sample>::iterator d_id;
                            if( reference_dependent )
                                {
                                    et_seq_search<seq_map,true> library_searcher( lib_idx, reference_counts, num_samples );

                                    auto seq_match = library_searcher.find( reads[ read_index ],
                                                                            std::get<2>( d_opts->seq_data ),
                                                                            seq_start,
                                                                            seq_length
                                                                          );
                                    
                                    //seq_duplicates[ seq_match->first.name ] = 0;

                                    if( seq_match != reference_counts.end() )
                                        {
                                            
                                            if( flexible_idx_data.size() > 1 )
                                                {
                                                    std::string concat_idx = "";
                                                    for( std::size_t curr_index = 0; curr_index < idx_match_list.size(); curr_index++ )
                                                        {
                                                            if( idx_match_list[curr_index] != index_map.end() )
                                                                {
                                                                    concat_idx.append( idx_match_list[curr_index]->first.seq );
                                                                }
                                                        }

                                                     d_id = index_map.find( sequence( "", concat_idx ) );
                                                }
                                            else
                                                {
                                                    d_id = index_map.find( sequence( "", idx_match_list[0]->first.seq ) );
                                                }

                                            if( d_id != index_map.end() )
                                                {
                                                    sample_id = d_id->second.id;
                                                    seq_match->second->at( sample_id ) += 1;
                                                    ++processed_success;
                                                    
                                                    if( !d_opts->diagnostic_fname.empty() )
                                                        {
                                                            
                                                            std::get<2>( diagnostic_map.find( d_id->second.name )->second[idx_match_list.size()] ) += 1;
                                                        }
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

                                    et_seq_search<seq_map,false> library_searcher( lib_idx, reference_counts, num_samples );


                                    if( flexible_idx_data.size() == 1 )
                                        {
                                            auto seq_match = library_searcher.find( reads[ read_index ],
                                                                                    std::get<2>( d_opts->seq_data ),
                                                                                    seq_start,
                                                                                    seq_length
                                                                                    );
                                            // check for duplicate sequence
                                            std::string seq_str = seq_match->first.name;
                                            /*if( seq_duplicates.find( seq_str ) != seq_duplicates.end() )
                                                {
                                                    seq_duplicates[seq_str]++;
                                                }
                                            else
                                                {
                                                    seq_duplicates.insert( std::make_pair( seq_str, 1 ) );
                                                }
                                            if( seq_match != reference_counts.end() )
                                                {
                                                    seq_match->second->at(0) += 1;
                                                    ++processed_success;
                                                }*/
                                            //seq_duplicates[ seq_str ] = 0;
                                        }
                                    else if( flexible_idx_data.size() > 1 )
                                        {
                                            std::string concat_idx = "";
                                            for( const auto& index : idx_match_list )
                                                {
                                                    sequential_map<sequence, sample> find_end;
                                                    if( index != find_end.end() )
                                                        {
                                                            concat_idx.append( index->first.seq );
                                                        }
                                                }

                                            d_id = index_map.find( sequence( "", concat_idx ) );

                                            if( d_id != index_map.end() )
                                                {
                                                    auto seq_match = library_searcher.find( reads[ read_index ],
                                                                                            std::get<2>( d_opts->seq_data ),
                                                                                            seq_start,
                                                                                            seq_length
                                                                                            );
                                                    // check for duplicate sequence
                                                    std::string seq_str = seq_match->first.name;
                                                    /*if( seq_duplicates.find( seq_str ) != seq_duplicates.end() )
                                                        {
                                                            seq_duplicates[seq_str]++;
                                                        }
                                                    else
                                                        {
                                                            seq_duplicates.insert( std::make_pair( seq_str, 1 ) );
                                                        }
                                                    if( seq_match != reference_counts.end() )
                                                        {
                                                            sample_id = d_id->second.id;
                                                            seq_match->second->at( sample_id ) += 1;
                                                            ++processed_success;
                                                        }
                                                        */

                                                    //seq_duplicates[seq_str] = 0;
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
                    idx_match_list.clear();
                }

            reads.clear();
            r2_seqs.clear();
            
        }
    total_time.stop();

    std::cout << processed_success << " records were found to be a match out of "
              << processed_total << " (" << ( (long double) processed_success / (long double) processed_total ) * 100
              << "%) successful.\n";
    // loop through index match totals
    for( std::size_t fif_index = 0; fif_index < flexible_idx_data.size(); ++fif_index )
        {
            std::cout << std::fixed << std::setprecision( 2 )
                    << ( ( long double ) index_match_totals[fif_index].second / ( long double ) processed_total ) * 100.00 << "% " <<
                                index_match_totals[fif_index].second << " of total records matched '" << flexible_idx_data[fif_index].idx_name << "'.\n";
        }
    if( reference_dependent )
        {
            std::cout << std::fixed << std::setprecision( 2 ) << ( ( long double ) processed_success / ( long double ) processed_total ) * 100.00 << "% " <<
                         processed_success << " of total records matched provided indexes + DNA tag.\n";
        }
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
    write_outputs( d_opts->output_fname, reference_counts, samplelist);
    write_diagnostic_output( d_opts, diagnostic_map);
}
                           //seq_duplicates

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

void module_demux::write_diagnostic_output( options_demux* d_opts, std::unordered_map<std::string,
             std::vector<std::tuple<std::string, std::string,std::size_t>>>& diagnostic_map )
{
    std::ofstream outfile( d_opts->diagnostic_fname, std::ios::out );

    std::string header = "Sample name\t";
    std::size_t index_right_offset = 1; // offset used to create a window allowing a concatenation to occur over each iteration.
    std::size_t num_indexes = diagnostic_map.begin()->second.size();
    // create column header for diagnostic file
    for( std::size_t curr_index = 0; curr_index < num_indexes; curr_index++ )
        {
            std::string header_col = "";
            for( std::size_t index_left_offset = 0; index_left_offset < index_right_offset ; index_left_offset++ )
                {
                    if( index_left_offset != 0 )
                        {
                            header_col.append( " + " );
                        }
                    std::string index_name = std::get<0>( diagnostic_map.begin()->second[index_left_offset] );
                    header_col.append(index_name);
                }
            header.append( header_col );
            header.append( "\t");
            ++index_right_offset;
        }
    header.append( "\n" );
    outfile << header;
    // output counts for index matches; sample major
    std::unordered_map<std::string,
             std::vector<std::tuple<std::string, std::string,std::size_t>>>::iterator sample_iter = diagnostic_map.begin();
    for( std::size_t curr_sample = 0; curr_sample < diagnostic_map.size(); curr_sample++ )
        {
            std::string line = "";
            line.append( sample_iter->first );
            for( const auto& index : sample_iter->second )
                {
                    line.append( "\t" );
                    line.append( std::to_string( std::get<2>( index ) ) );
                }
            outfile << line;
            outfile << "\n";
            ++sample_iter;
        }

    outfile << "\n";
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

    std::map<std::string, std::size_t> dups;

    for( index = 0; index < seq_scores.size(); ++index )
        {
            const sequence& curr = seq_iter->first;
            const std::vector<std::size_t> *curr_counts = seq_iter->second;
                    if( ( ++dups [ curr.seq ] ) == 1)
                        {
                            outfile << curr.name << DELIMITER;
                

                            for( second_index = 0; second_index < curr_counts->size() - 1; ++second_index )
                                {
                                    outfile << curr_counts->at( second_index ) << DELIMITER;
                                }
                            outfile << curr_counts->at( curr_counts->size() - 1 ) << NEWLINE;
                        }
            ++seq_iter;
            delete curr_counts;
        }
    outfile.close();

    for( const auto row : dups )
        {
            if( row.second > 1 )
                {
                    std::cout << "Duplicate: " << row.first << ", " << row.second;
                }
        }
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
                                     std::vector<sample>& samplelist,
                                     sequential_map<sequence, sample>& seq_lookup
                                   )
{
    std::unordered_map<std::string, std::string> updated_index_seqs;    
    for( std::size_t index = 0; index < index_seqs.size(); ++index )
        {
            bool index_seq_found = false;
            std::size_t curr_sample = 0;
            // while the current index seq id has not been found in a sample and while there are still samples
            while( !index_seq_found && curr_sample < samplelist.size() )
                {
                    // look at each string id and compare
                    for( std::size_t curr_id = 0; curr_id < samplelist[curr_sample].string_ids.size(); ++curr_id )
                        {
                            // if the current id does match the index seqs id
                            if( index_seqs[index].name.compare( samplelist[curr_sample].string_ids[curr_id] ) == 0 )
                                {
                                    // mark as the index seq was found
                                    index_seq_found = true;
                                    // add the index seq to the updated list (order doesnt matter)
                                    updated_index_seqs.emplace( index_seqs[index].name, index_seqs[index].seq );
                                    seq_lookup[sequence( "", index_seqs[index].seq )] = samplelist[curr_sample];
                                }
                        }
                    curr_sample++;
                }
        }
    // future function: add_match_references -> add the partial matches and complete matches by creating 
    for( std::size_t curr_sample = 0; curr_sample < samplelist.size(); ++curr_sample )
        {
            std::string concat_seq = "";
            for( std::size_t curr_idx = 0; curr_idx < samplelist[curr_sample].string_ids.size(); ++curr_idx )
                {
                    // obtain the barcode for the current index for the current sample and add to sequence concat to create element.
                    auto sample_seq = updated_index_seqs.find( samplelist[curr_sample].string_ids[curr_idx] )->second;
                    concat_seq.append( sample_seq );
                    map[ sequence( "", concat_seq ) ] = samplelist[ curr_sample ];
                }
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

