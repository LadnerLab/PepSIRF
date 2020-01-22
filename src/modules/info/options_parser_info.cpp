#include "options_parser_info.h"
#include <boost/program_options.hpp>

options_parser_info::options_parser_info() = default;
options_parser_info::~options_parser_info() = default;

bool options_parser_info
::parse( int argc, char ***argv, options *opts )
{
    options_info *opts_info = (options_info*) opts;
    namespace po = boost::program_options;
    po::variables_map vm;

    po::options_description desc( "PepSIRF: Peptide-based Serological Immune "
                                  "Response Framework info module.\n"
                                );

    desc.add_options()
        (
         "help,h", "Produce help message and exit.\n"
         "This module is used to gather information about a score matrix. "
         "By default, the number of samples and peptides in the matrix will be output. "
         "Additional flags may be used to output different information.\n"
        )



        ;

    po::store( po::command_line_parser( argc, *argv ).options( desc ).run(), vm);

    if( vm.count( "help" ) 
        || argc == 2 
        )
        {
            std::cout << desc << std::endl;
            return false;
        }
    else
        {
            po::notify( vm );
            return true;
        }

    return true;
}
