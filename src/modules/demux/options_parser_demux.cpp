#include "options_parser_demux.h"

options_parser_demux::options_parser_demux() = default;

void options_parser_demux::parse( int argc, char ***argv, options *opts )
{
    options_demux *opts_demux = (options_demux*) opts;
    namespace po = boost::program_options;
    po::variables_map vm;

    po::options_description desc( "PepSIRF: Peptide-based Serological Immune Response Framework" );
    desc.add_options()
        ( "help,h", "Produce help message" )
        ( "input_r1", po::value<std::string>( &opts_demux->input_r1_fname )->required(), "Input forward reads fastq file to parse.")
        ( "library,l", po::value<std::string>( &opts_demux->library_fname )->required(), "Designed library containing nucleic acid peptides. "
                                                             "Library should be in fasta format and should contain "
                                                             "sequences that were used to design input_r1."
        )
        ( "read_per_loop,r", po::value<long int>( &opts_demux->read_per_loop )->default_value( opts_demux->DEFAULT_READ_PER_LOOP, "The number of fastq "
          "records read a time. A higher value will result in more memory usage by the program, but will also result in fewer disk accesses, "
          "increasing performance of the program." )
        )
        ( "num_threads,t", po::value<int>( &opts_demux->num_threads )->default_value( opts_demux->DEFAULT_NUM_THREADS ), "Number of threads to use for analyses." );

    po::store( po::command_line_parser( argc, *argv ).options( desc ).allow_unregistered().run(), vm);

    if( vm.count( "help" ) )
        {
            std::cout << desc << std::endl;
        }
    else
        {
            po::notify( vm );
        }
}
