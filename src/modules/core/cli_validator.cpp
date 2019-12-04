#include "cli_validator.h"
#include <iostream>
#include <algorithm>
#include <string>
#include <unordered_set>
#include "modules.h"

cli_validator::cli_validator() = default;

bool cli_validator::validate( int argc, char ***argv )
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
                    std::cout << "The currently available modules are:\n";
                    std::cout << modules::join_names( " - ", "\n" );
                    std::cout << "\n";
                    std::cout << "--help, -h displays this message, while "
                                 "'pep_sirf module_name --help' will display "
                                 "the help for " 
                                 "the module module_name.\n";
                }
            return false;
        }

    if( modules::is_module( argv_loc[ 2 ] ) )
        {
            return true;
        }
    throw std::runtime_error( "Invalid module name entered" );
}
