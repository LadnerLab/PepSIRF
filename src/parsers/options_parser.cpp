#include "options_parser.h"

OptionsParser::OptionsParser() = default;

void OptionsParser::parse( int argc, char ***argv )
{
    namespace po = boost::program_options;
    char **argv_loc = *argv;

    try
        {
            po::options_description desc( "Allowed options" );
            desc.add_options()
                ( "help", "Produce help message" );

            po::variables_map vm;
            po::store( po::parse_command_line( argc, argv_loc, desc ), vm );

            if( vm.count( "help" ) )
                {
                    std::cout << desc << std::endl;
                }

        }
    catch( std::exception& e )
        {
            std::cerr << "Error: " << e.what() << std::endl;
        }
    catch( ... )
        {
            std::cerr << "Unknown exception encountered!" << std::endl;
        }
}
