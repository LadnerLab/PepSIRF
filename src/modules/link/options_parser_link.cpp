#include "options_parser_link.h"
#include "options_link.h"

#include <boost/program_options.hpp>

bool options_parser_link::parse( int argc, char ***argv, options *opts )
{
    options_link *opts_bin = (options_link*) opts;
    namespace po = boost::program_options;
    po::variables_map vm;

    po::options_description desc( "PepSIRF: Peptide-based Serological Immune "
                                  "Response Framework link module.\n"
                                  );

    desc.add_options()
        (
         "help,h", "Produce help message and exit.\n"
        );

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
}
