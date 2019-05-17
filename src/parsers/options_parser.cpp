#include "options_parser.h"

options_parser::options_parser() = default;

options options_parser::parse( int argc, char ***argv )
{
    namespace po = boost::program_options;
    char **argv_loc = *argv;

    po::options_description desc( "Allowed options" );
    desc.add_options()
        ( "help", "Produce help message" );

    po::variables_map vm;
    po::store( po::parse_command_line( argc, argv_loc, desc ), vm );

    if( vm.count( "help" ) )
        {
            std::cout << desc << std::endl;
        }

    return options();
}
