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
                    std::cout << "\nUSAGE: pepsirf [ --help | module_name <module_args*> ] " << "\n";
                    std::cout << "The currently available modules are:\n";
                    std::cout << " -demux\n";
                    std::cout << " -norm\n";
                    std::cout << " -deconv\n\n";
                    std::cout << "--help, -h displays this message, while 'pepsirf module_name --help' will display the help for " 
                        "the module module_name.\n";
                    return nullptr;
                }

            if( !arg.compare( "demux" ) )
                {
                    return new options_parser_demux();
                }
            else if( !arg.compare( "deconv" ) )
                {
                    return new options_parser_deconv();
                }
            else if( !arg.compare( "norm" ) )
                {
                    return new options_parser_normalize();
                }

        }
    throw std::runtime_error( "Invalid module name entered" );
}
