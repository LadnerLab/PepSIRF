#include "module_demux.h"

module_demux::module_demux() 
{
    name = "Demux";
}

void module_demux::run( options *opts )
{
    // options from the command line
    options_demux *d_opts = (options_demux*) opts;

    struct time_keep::timer total_time;
    std::unordered_map<sequence, std::vector<std::size_t>> reference_counts;

    omp_set_num_threads( opts->num_threads );

    total_time.start = omp_get_wtime();

    // vector to store the .fna sequences that represent a designed library
    std::vector<sequence> library_seqs;
    std::vector<sequence> reads;

    reads.reserve( d_opts->read_per_loop );

    // create parsers for the encoded library and the fastq reads.
    fasta_parser fasta_p;
    fastq_parser fastq_p;

    // open the read file so we can iterate over it 
    std::ifstream reads_file( d_opts->input_r1_fname, std::ios_base::in );

    library_seqs = fasta_p.parse( d_opts->library_fname );
    add_seqs_to_map( reference_counts, library_seqs, 11 );

    // while we still have reads to process
    while( fastq_p.parse( reads_file, reads, d_opts->read_per_loop  ) )
        {
            reads.clear();
        }

    total_time.end = omp_get_wtime();

    std::cout << "Total elapsed time (seconds): " << time_keep::get_elapsed( total_time ) << ".\n";
}


std::string module_demux::get_name()
{
    return name;
}

void module_demux::add_seqs_to_map( std::unordered_map<sequence, std::vector<std::size_t>>& input_map, std::vector<sequence>& seqs, size_t num_samples )
{
    std::size_t index = 0;
    input_map.reserve( seqs.size() );

    for( index = 0; index < seqs.size(); ++index )
        {
            input_map[ seqs[ index ] ] = std::vector<std::size_t>( num_samples );
        }
}
