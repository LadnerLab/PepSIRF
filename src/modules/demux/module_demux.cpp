#include <sstream>
#include <stdexcept>

#include "file_io.h"
#include "module_demux.h"
#include "fastq_sequence.h"
#include "fastq_score.h"
#include "et_search.h"
#include "nt_aa_translator.h"
#include <thread>
#include <mutex>
#include "fs_tools.h"

module_demux::module_demux()
{
    name = "Demux";
}

void module_demux::run( options *opts )
{
    // options from the command line
    options_demux *d_opts = (options_demux*) opts;

    std::string index_str;
    std::mutex mtx;

    std::map<std::string, std::size_t> duplicate_map;

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
            if( !d_opts->input_r2_fname.empty() )
                {
                    flexible_idx_data.emplace_back( d_opts->sample_indexes[1],
                                                    "r2",
                                                    std::get<0>( d_opts->index2_data ),
                                                    std::get<1>( d_opts->index2_data ),
                                                    std::get<2>( d_opts->index2_data ));
                }
        }
    else
        {
            fif_parser flex_parser;
            flexible_idx_data = flex_parser.parse( d_opts->flexible_idx_fname );
            d_opts->sample_indexes.clear();
            for( const auto& flex_index : flexible_idx_data )
                {
                    d_opts->sample_indexes.emplace_back(flex_index.idx_name);
                }
        }


    std::size_t read_index = 0;
    struct time_keep::timer total_time;
    parallel_map<sequence, std::vector<std::size_t>*> reference_counts;
    parallel_map<sequence, std::size_t> non_perfect_match_seqs;
    // std::unordered_map<std::string, std::vector<std::tuple<std::string, std::string, std::size_t>>> diagnostic_map;

    omp_set_num_threads( opts->num_threads );

    total_time.start();

    // vector to store the .fna sequences that represent a designed library
    std::vector<sequence> library_seqs;
    std::vector<sequence> dna_tags;
    std::vector<fastq_sequence> r2_seqs;
    std::vector<fastq_sequence> reads;

    reads.reserve( d_opts->read_per_loop );
    // create parsers for the encoded library and the fastq reads.
    fasta_parser fasta_p;
    fastq_parser fastq_p;
    // parse samplelist
    samplelist_parser samplelist_p;
    std::vector<sample> samplelist = samplelist_p.parse( d_opts, flexible_idx_data );

    std::ifstream reads_file( d_opts->input_r1_fname, std::ios_base::in );
    std::ifstream r2_reads;

    bool r2_reads_included = d_opts->input_r2_fname.length() > 0;

    if( r2_reads_included )
        {
            r2_reads.open( d_opts->input_r2_fname, std::ios_base::in );
        }
    else
        {
            // Should we give a warning for the case of how many index col names are given?
            if( d_opts->sample_indexes.size() > 1 )
                {
                    std::cout << "WARNING: \"--input_r2\" does not include a filepath while more "
                    "than one samplelist index column name has been given. See the \"--sname\" "
                    "and \"--samplelist\" options for information on default values and usage.\n";
                }
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
    sequential_map<sequence, sample> index_map;
    sequential_map<sequence, sample> seq_lookup;
    dna_tags = fasta_p.parse( d_opts->index_fname );
    create_index_map( index_map, dna_tags, samplelist, seq_lookup);
    std::size_t count = 0;
    for( const auto &index_seq : index_map )
    {
        auto seq = index_seq.first;
        auto sample = index_seq.second;
        count++;
    }
    
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
    index_idx.index( dna_tags );

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
    phmap::parallel_flat_hash_map<sample, std::vector<std::size_t>> diagnostic_map;
    if( !d_opts->diagnostic_fname.empty() )
        {
            create_diagnostic_map( reference_dependent, diagnostic_map, samplelist);
        }

    if( !d_opts->fastq_out.empty() )
        {
            auto output_path = fs_tools::path( d_opts->fastq_out );
            bool dir_exists = !fs_tools::create_directories( output_path );

            if( dir_exists )
            {
                std::cout << "WARNING: the directory '" << d_opts->fastq_out
                          << "' exists, any files with "
                          << "colliding filenames will be overwritten!\n";
            }

            std::stringstream build;

            for( auto sample : samplelist ) 
                {
                    std::ofstream clear_file;
                    build << d_opts->fastq_out << "/" << sample.name << ".fastq";
                    clear_file.open(build.str(), std::ios_base::out);
                    clear_file.close();
                    build.clear();
                }
        }


    while( fastq_p.parse( reads_file_ref, reads, d_opts->read_per_loop  ) )
        {

            if( r2_reads_included )
                {
                    fastq_p.parse( r2_reads_ref, r2_seqs, d_opts->read_per_loop );
                }
            std::map<std::string, std::vector<fastq_sequence>> fastq_output;

#ifndef __clang__

           #pragma omp parallel for private( seq_iter, nuc_seq, read_index, index_str, adapter, sample_id, duplicate_map,  \
                                              idx_match_list) \
               shared( seq_start, seq_length, d_opts, num_samples, reference_counts, library_seqs, dna_tags, r2_seqs) \
               reduction( +:processed_total, processed_success, concatemer_found ) \


            for( read_index = 0; read_index < reads.size(); ++read_index )
#endif //__clang__
#ifdef __clang__
            for( read_index = 0; read_index < reads.size(); ++read_index )
#endif
                {
                    // Identify found matches
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
#ifndef __clang__
                                                    #pragma omp critical
                                                        {
#endif
                                                            index_match_totals[curr_index].second++;
                                                            if( !d_opts->diagnostic_fname.empty() )
                                                            {
                                                                // for each sample/index count element
                                                                for( auto& sample_index_count : diagnostic_map )
                                                                    {
                                                                        // while matches are found, increment to current index
                                                                        std::size_t review_index = 0;
                                                                        bool has_match = true;
                                                                        // for each found match,
                                                                        // verify up to current index that all matches are valid
                                                                        while( has_match && review_index < curr_index )
                                                                            {
                                                                                if( sample_index_count.first.string_ids[review_index].compare(matched_ids[review_index]) != 0 )
                                                                                    {
                                                                                        has_match = !has_match;
                                                                                    }
                                                                                review_index++;
                                                                            }
                                                                        // increment count if indexes up to current match
                                                                        if( has_match
                                                                            && sample_index_count.first.string_ids[review_index].compare(matched_ids[review_index]) == 0 )
                                                                            {
                                                                                sample_index_count.second[review_index] += 1;
                                                                            }
                                                                    }
                                                            }
#ifndef __clang__
                                                        }
#endif
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
                    // Add valid matches
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
                    // Handle match found
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
                                    if( seq_match != reference_counts.end() )
                                        {
                                            if( flexible_idx_data.size() > 1 )
                                                {
                                                    std::string concat_idx = "";

#ifndef __clang__
                                                    #pragma omp critical
                                                    {
#endif
                                                    for( std::size_t curr_index = 0; curr_index < idx_match_list.size(); curr_index++ )
                                                        {
                                                            if( idx_match_list[curr_index] != index_map.end() )
                                                                {
                                                                    concat_idx.append( idx_match_list[curr_index]->first.seq );
                                                                }
                                                        }

                                                     d_id = index_map.find( sequence( "", concat_idx ) );
#ifndef __clang__
                                                    }
#endif
                                                }
                                            else
                                                {
                                                    d_id = index_map.find( sequence( "", idx_match_list[0]->first.seq ) );
                                                }

                                            if( d_id != index_map.end() )
                                                {
                                                    sample_id = d_id->second.id;

                                                    if( !d_opts->fastq_out.empty() )
                                                        {
#ifndef __clang__
                                                            #pragma omp critical
                                                            {
#endif                                                          
                                                                fastq_output[d_id->second.name].emplace_back(reads[ read_index ]);
#ifndef __clang__
                                                            }
#endif
                                                        }
#ifndef __clang__
                                                    #pragma omp critical
                                                        {
#endif
                                                            seq_match->second->at( sample_id ) += 1;
#ifndef __clang__
                                                        }
#endif
                                                    ++processed_success;
                                                    
                                                    if( !d_opts->diagnostic_fname.empty() )
                                                        {
#ifndef __clang__
                                                            #pragma omp critical
                                                                {
#endif
                                                                    diagnostic_map.find( d_id->second )->second[idx_match_list.size()] += 1;
#ifndef __clang__
                                                                }
#endif
                                                            
                                                        }
                                                }
                                        }
                                    else if( seq_match == reference_counts.end()
                                             && found_concatemer() )
                                        {
#ifndef __clang__
                                            #pragma omp critical
                                                {
#endif
                                                    ++concatemer_found;
#ifndef __clang__
                                                }
#endif
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
                                            
                                            
                                            if( seq_match != reference_counts.end() )
                                                {
                                                    // if seq_match found, increase count for given sample
                                                    auto sample_id = idx_match_list[0]->second.id;

#ifndef __clang__
                                                    #pragma omp critical
                                                    {
#endif
                                                        fastq_output[d_id->second.name].emplace_back(reads[ read_index ]);
#ifndef __clang__
                                                    }
#endif
#ifndef __clang__
                                                    #pragma omp critical
                                                        {
#endif
                                                            seq_match->second->at(sample_id) += 1;
#ifndef __clang__
                                                        }
#endif
                                                    ++processed_success;
                                                }        
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
                                                    if( seq_match != reference_counts.end() )
                                                        {
                                                            sample_id = d_id->second.id;
#ifndef __clang__
                                                            #pragma omp critical
                                                                {
#endif
                                                                    seq_match->second->at( sample_id ) += 1;
                                                                    //ISS169
#ifndef __clang__
                                                                }
#endif

                                                            ++processed_success;
                                                        }
                                                }
                                        }
                                    else if( found_concatemer() )
                                        {
#ifndef __clang__
                                            #pragma omp critical
                                                {
#endif
                                                    ++concatemer_found;
#ifndef __clang__
                                                }
#endif
                                        }

                                }

                        }
                    // record the number of records that are processed
                    ++processed_total;
                    idx_match_list.clear();
                }

#ifndef __clang__
            #pragma omp critical
                {
#endif
                    write_fastq_output(fastq_output, d_opts->fastq_out);
#ifndef __clang__
                }
#endif

            fastq_output.clear();

            reads.clear();
            r2_seqs.clear();
        }
    total_time.stop();
    // check for duplicates
#ifndef __clang__
    #pragma omp critical
    {
#endif
    for( const auto& library_seq : library_seqs )
    {
        ++duplicate_map[library_seq.seq];
    }
#ifndef __clang__
    }
#endif

    std::cout << processed_success << " records were found to be a match out of "
              << processed_total << " (" << ( (long double) processed_success / (long double) processed_total ) * 100
              << "%) successful.\n";
    // loop through index match totals
    std::string match_reference = "";
    for( std::size_t fif_index = 0; fif_index < flexible_idx_data.size(); ++fif_index )
        {
            match_reference.append( flexible_idx_data[fif_index].idx_name );
            std::cout << std::fixed << std::setprecision( 2 )
                    << ( ( long double ) index_match_totals[fif_index].second / ( long double ) processed_total ) * 100.00 << "% " <<
                                index_match_totals[fif_index].second << " of total records matched '" << match_reference << "'.\n";
            match_reference.append( " + " );
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

    if( !d_opts->diagnostic_fname.empty() )
        {
            write_diagnostic_output( d_opts, diagnostic_map);
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
                           duplicate_map,
                           !d_opts->library_fname.empty(),
                           samplelist
                         );
        }
    write_outputs( d_opts->output_fname, reference_counts, duplicate_map, !d_opts->library_fname.empty(), samplelist);
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

            std::string copy_name = strs[ 0 ];
            std::string copy_seq = iter->first.seq;
            current.set_name(copy_name);
            current.set_seq(copy_seq);

            //std::cout << current.name << std::endl;

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
                    // count occurrences of distinct sequences
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

    #pragma omp parallel for private( index ) shared ( seqs, input_map, num_samples )
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

void module_demux::create_diagnostic_map( bool reference_dependent,
                                          phmap::parallel_flat_hash_map<sample,std::vector<std::size_t>>& diagnostic_map,
                                          std::vector<sample> samplelist )
{
    std::size_t column_size = reference_dependent ? samplelist[0].string_ids.size() + 1 : samplelist[0].string_ids.size();
    for( const auto &sample : samplelist )
        {
            diagnostic_map.emplace( sample, std::vector<std::size_t>(column_size, 0) );
        }

}

void module_demux::write_diagnostic_output( options_demux* d_opts, phmap::parallel_flat_hash_map<sample, std::vector<std::size_t>>& diagnostic_map )
{
    std::ofstream outfile( d_opts->diagnostic_fname, std::ios::out );
    outfile << "Sample name\t";
    std::string header_column = "";
    
    // create column header for diagnostic file
    // for each column name
    for( std::size_t curr_column = 0; curr_column < d_opts->sample_indexes.size(); curr_column++ )
        {
            std::size_t last_index_name = curr_column + 1;
            std::size_t curr_index_name = 0;
            // for each index name - with index end at column + 1
            while( curr_index_name < last_index_name )
                {
                    // add index name
                    header_column += d_opts->sample_indexes[curr_index_name];
                    if( curr_index_name + 1 < last_index_name )
                        {
                            header_column += " + ";
                        }
                    curr_index_name++;
                }
            header_column += "\t";
        }

    outfile << header_column;
    if( !d_opts->library_fname.empty() )
        {
            outfile << "DNA tags";
        }
    outfile << "\n";
    // output counts for index matches; sample major
    phmap::parallel_flat_hash_map<sample, std::vector<std::size_t>>::iterator sample_iter = diagnostic_map.begin();
    for( std::size_t curr_sample = 0; curr_sample < diagnostic_map.size(); curr_sample++ )
        {
            std::string line = "";
            line += sample_iter->first.name;
            for( std::size_t curr_index = 0; curr_index < sample_iter->second.size(); curr_index++ )
                {
                    line += "\t" + std::to_string(sample_iter->second[curr_index]);
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
                                  std::map<std::string, std::size_t> duplicate_map,
                                  bool ref_dependent,
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

            if(!ref_dependent || duplicate_map[curr.seq] == 1)
                {
                     outfile << curr.name << DELIMITER;

                    for( second_index = 0; second_index < samples.size() - 1; ++second_index )
                        {
                            outfile << curr_counts->at( second_index ) << DELIMITER;
                        }
                    outfile << curr_counts->at( samples.size() - 1 ) << NEWLINE;
                }

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

void module_demux::create_index_map( sequential_map<sequence, sample>& index_map,
                                     std::vector<sequence>& dna_tags,
                                     std::vector<sample>& samplelist,
                                     sequential_map<sequence, sample>& seq_lookup
                                   )
{
    for( const auto& sample : samplelist )
    {
        std::string concat_sequence = "";
        for( size_t sample_id_index = 0; sample_id_index < sample.string_ids.size(); sample_id_index++ )
        {
            for( size_t dna_tag_index = 0; dna_tag_index < dna_tags.size(); dna_tag_index++ )
            {
                if( dna_tags[dna_tag_index].name.compare(sample.string_ids[sample_id_index]) == 0 )
                {
                    seq_lookup.emplace( sequence( "", dna_tags[dna_tag_index].seq ), sample );
                    concat_sequence += dna_tags[dna_tag_index].seq;
                    index_map.emplace( sequence( "", concat_sequence ), sample );
                    break;
                }
            }
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