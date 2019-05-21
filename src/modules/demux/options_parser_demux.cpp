#include "options_parser_demux.h"

options_parser_demux::options_parser_demux() = default;

void options_parser_demux::parse( int argc, char ***argv, options_demux& opts )
{
    namespace po = boost::program_options;
    char **argv_loc = *argv;
    po::variables_map vm;

    po::options_description desc( "PepSIRF: Peptide-based Serological Immune Response Framework" );
    desc.add_options()
        ( "help,h", "Produce help message" )
        ( "input_r1", po::value<std::string>( &opts.input_r1_fname )->required(), "Input forward reads fastq file to parse.")
        ( "input_r2", po::value<std::string>( &opts.input_r2_fname )->required(), "Input reverse reads fastq file to parse.")
        ( "library,l", po::value<std::string>( &opts.library_fname )->required(), "Designed library containing amino acid peptides. "
                                                             "Library should be in fasta format and should contain "
                                                             "sequences that were used to design input_r1 and input_r2."
        )
        ( "num_threads,t", po::value<int>( &opts.num_threads )->default_value( opts.DEFAULT_NUM_THREADS ), "Number of threads to use for analyses." );

    po::store( po::parse_command_line( argc, argv_loc, desc ), vm );

    if( vm.count( "help" ) )
        {
            std::cout << desc << std::endl;
        }

    po::notify( vm );
}
