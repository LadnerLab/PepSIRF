#include "options_parser_factory.h"

options_parser_factory::options_parser_factory() = default;

options_parser *options_parser_factory::create( int argc, char ***argv )
{
    namespace po = boost::program_options;
    char **argv_loc = *argv;
    std::string mod_name;
    po::variables_map vm;

    po::options_description desc( "PepSIRF: Peptide-based Serological Immune Response Framework" );
    desc.add_options()
        ( "help,h", "Produce help message" )
        ( "mod_name", po::value<std::string>(), "Name of module to run.");


    po::positional_options_description pos;
    po::options_description all_options;

    all_options.add( desc );

    pos.add( "mod_name", 1 );
    po::store( po::command_line_parser( argc, argv_loc ).
               options( all_options ).
               positional( pos ).
               allow_unregistered().
               run() , vm
             );


    mod_name = vm[ "mod_name" ].as<std::string>();

    std::transform( mod_name.begin(), mod_name.end(), mod_name.begin(), ::tolower );

    if( mod_name.compare( "demux" ) == 0 )
        {
            return new options_parser_demux();
        }
    throw std::runtime_error( "Invalid module name entered" );
}
