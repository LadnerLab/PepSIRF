#include <boost/program_options.hpp>
#include <iostream>
#include <iterator>
#include <exception>
#include <stdexcept>

#include "options_parser_deconv.h"


options_parser_deconv::options_parser_deconv() = default;

bool options_parser_deconv::parse( int argc, char ***argv, options *opts )
{
    options_deconv *opts_deconv = (options_deconv*) opts;
    namespace po = boost::program_options;
    po::variables_map vm;

    po::options_description desc( "PepSIRF: Peptide-based Serological Immune Response Framework species deconvolution module. \n"
                                );
    desc.add_options()
        ( "help,h", "Produce help message" );

    po::store( po::command_line_parser( argc, *argv ).options( desc ).allow_unregistered().run(), vm);

    if( vm.count( "help" ) )
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
