#include "options_parser.h"

options_parser::options_parser() = default;
options_parser::~options_parser() = default;

void options_parser::parse( int argc, char ***argv, options *opts )
{
    namespace po = boost::program_options;
    char **argv_loc = *argv;
    po::variables_map vm;

    po::options_description desc( "PepSIRF: Peptide-based Serological Immune Response Framework" );
    desc.add_options()
        ( "help,h", "Produce help message" )
        ( "num_threads,t", po::value<int>( &opts->num_threads )->default_value( opts->DEFAULT_NUM_THREADS ), "Number of threads to use for analyses." );

    po::store( po::parse_command_line( argc, argv_loc, desc ), vm );

    if( vm.count( "help" ) )
        {
            std::cout << desc << std::endl;
        }

    po::notify( vm );
}
