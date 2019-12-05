#include "options_parser_subjoin.h"
options_parser_subjoin::options_parser_subjoin() = default;


bool options_parser_subjoin::parse( int argc, char ***argv, options *opts )
{
    options_subjoin *opts_subjoin = (options_subjoin*) opts;
    namespace po = boost::program_options;
    po::variables_map vm;


    po::options_description desc( "PepSIRF: Peptide-based Serological Immune "
                                  "Response Framework subjoin module"
                                );

    desc.add_options()
        (
         "help,h", "Produce help message\n"
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
