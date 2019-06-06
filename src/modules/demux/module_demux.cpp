#include "module_demux.h"


module_demux::module_demux() 
{
    name = "Demux";
}

void module_demux::run( options *opts )
{
    // options from the command line
    options_demux *d_opts = (options_demux*) opts;

    std::string index_str;

    int forward_start  = 12;
    int forward_length = 12;

    int reverse_start    = 21;
    int reverse_length   = 23;

    size_t read_index = 0;

    struct time_keep::timer total_time;
    parallel_map<sequence, std::vector<std::size_t>*> reference_counts;
    sequential_map<std::string, sample> index_map;

    omp_set_num_threads( opts->num_threads );

    total_time.start = omp_get_wtime();

    // vector to store the .fna sequences that represent a designed library
    std::vector<sequence> library_seqs;
    std::vector<sequence> index_seqs;
    std::vector<sequence> reads;

    reads.reserve( d_opts->read_per_loop );

    // create parsers for the encoded library and the fastq reads.
    fasta_parser fasta_p;
    fastq_parser fastq_p;

    samplelist_parser samplelist_p;
    std::vector<sample> samplelist = samplelist_p.parse( d_opts->samplelist_fname );

    // open the read file so we can iterate over it 
    std::ifstream reads_file( d_opts->input_r1_fname, std::ios_base::in );
    std::ifstream r2_reads;

    if( d_opts->input_r2_fname.length() > 0 )
        {
            r2_reads.open( d_opts->input_r2_fname, std::ios_base::in );
        }

    library_seqs = fasta_p.parse( d_opts->library_fname );
    index_seqs   = fasta_p.parse( d_opts->index_fname );

    demux::create_index_map( index_map, index_seqs, samplelist );
    sequential_map<std::string, sample> sample_table;

    add_seqs_to_map( reference_counts, library_seqs, samplelist.size() );

    std::string adapter;
    std::string nuc_seq;

    parallel_map<sequence, std::vector<std::size_t>*>::iterator seq_iter;

    unsigned int seq_start = 43;
    unsigned int seq_length = 90;

    std::size_t processed_total   = 0;
    std::size_t processed_success = 0;

    // while we still have reads to process, loop terminates when no fastq records
    // are read from the file
    while( fastq_p.parse( reads_file, reads, d_opts->read_per_loop  ) )
        {

            if( d_opts->input_r2_fname.length() > 0 )
                {
                    // get the index sequences for this set of sequences
                    fastq_p.parse( r2_reads, index_seqs, d_opts->read_per_loop );
                }

            #pragma omp parallel for private( seq_iter, nuc_seq, read_index, index_str, adapter ) \
                shared( seq_start, seq_length, d_opts, reference_counts, library_seqs, index_seqs ) reduction( +:processed_total, processed_success )
            for( read_index = 0; read_index < reads.size(); ++read_index )
                {
                    // get the forward and the reverse indexes from the sequence, grab the id
                    // for this index. This gives the index of the location at the sequence to increment
                    adapter = reads[ read_index ].seq.substr( forward_start, forward_length );

                    if( index_map.find( adapter ) != index_map.end() )
                        {
                            nuc_seq = reads[ read_index ].seq.substr( seq_start, seq_length );
                            seq_iter = reference_counts.find( sequence( "", nuc_seq ) );

                            if( seq_iter != reference_counts.end() )
                                {
                                    ( seq_iter->second->at( index_map[ adapter ].id ) )++;

                                    ++processed_success;
                                }
                        }

                    // record the number of records that are processed
                    ++processed_total;
                }
            reads.clear();
        }

    total_time.end = omp_get_wtime();

    if( d_opts->input_r2_fname.length() > 0 )
        {
            r2_reads.close();
        }

    std::cout << "Processed " << processed_total << " records in " << time_keep::get_elapsed( total_time ) << " seconds.\n";
    std::cout << processed_success << " records were an exact match out of "
              << processed_total << " (" << ( (long double) processed_success / (long double) processed_total ) * 100
              << "%) successful.\n";

    write_outputs( "out_test", reference_counts, samplelist );

}


std::string module_demux::get_name()
{
    return name;
}

void module_demux::add_seqs_to_map( parallel_map<sequence, std::vector<std::size_t>*>& input_map, std::vector<sequence>& seqs, size_t num_samples )
{
    std::size_t index        = 0;
    std::size_t inner_index  = 0;

    input_map.reserve( seqs.size() );

    #pragma omp parallel for private( index, inner_index  ) shared ( seqs, input_map, num_samples )
    for( index = 0; index < seqs.size(); ++index )
        {
            input_map[ seqs[ index ] ] = new std::vector<std::size_t>( num_samples );

            for( inner_index = 0; inner_index < num_samples; ++inner_index )
                {
                    input_map[ seqs[ index ] ]->at( inner_index ) = 0;
                }
        }
}

void module_demux::write_outputs( std::string outfile_name,
                                  parallel_map<sequence, std::vector<std::size_t>*>& seq_scores,
                                  std::vector<sample>& samples
                                )
{
    std::ofstream outfile( outfile_name, std::ofstream::out );
    parallel_map<sequence, std::vector<std::size_t>*>::iterator seq_iter = seq_scores.begin();


    std::size_t index        = 0;
    std::size_t second_index = 0;

    for( index = 0; index < seq_scores.size(); ++index )
        {
            const sequence& curr = seq_iter->first;
            const std::vector<std::size_t> *curr_counts = seq_iter->second;

            outfile << curr.name << ",";

            for( second_index = 0; second_index < samples.size() - 1; ++second_index )
                {
                    outfile << curr_counts->at( second_index ) << ',';
                }
            outfile << curr_counts->at( samples.size() - 1 ) << "\n";

            ++seq_iter;
            delete curr_counts;
        }
    outfile.close();
}

void demux::create_index_map( sequential_map<std::string, sample>& map,
                              std::vector<sequence>& index_seqs,
                              std::vector<sample>& samplelist
                            )
{
    sequential_map<std::string, sample> sample_table;

    unsigned int index = 0;

    for( index = 0; index < samplelist.size(); ++index )
        {
            sample_table[ samplelist[ index ].get_ids() ] = samplelist[ index ];
        }

    for( index = 0; index < index_seqs.size(); ++index )
        {
            if( sample_table.find( index_seqs[ index ].name ) != sample_table.end() )
                {
                    map[ index_seqs[ index ].seq ] = sample_table[ index_seqs[ index ].name ];
                }
        }
}
