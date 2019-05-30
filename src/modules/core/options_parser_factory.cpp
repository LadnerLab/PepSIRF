#include "options_parser_factory.h"

options_parser_factory::options_parser_factory() = default;

options_parser *options_parser_factory::create( int argc, char ***argv )
{
    char **argv_loc = *argv;
 
    std::string desc = "PepSIRF: Peptide-based Serological Immune Response Framework";
    std::string arg;
    if( argc >= 2 )
        {
            arg = argv_loc[ 1 ];
            std::transform( arg.begin(), arg.end(), arg.begin(), ::tolower );

            if( !arg.compare( "-h" ) ||
                !arg.compare( "--help" )
              )
                {
                    std::cout << desc << "\n";
                    std::cout << "\nUSAGE: pep_sirf [ --help | module_name <module_args*> ] " << "\n";
                    std::cout << "--help, -h displays this message, while 'pep_sirf module_name --help' will display the help for " 
                        "the module module_name.\n";
                    return nullptr;
                }

            if( !arg.compare( "demux" ) )
                {
                    return new options_parser_demux();
                }
        }
    throw std::runtime_error( "Invalid module name entered" );
}
