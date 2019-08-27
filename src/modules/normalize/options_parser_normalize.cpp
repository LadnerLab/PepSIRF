#include <boost/program_options.hpp>
#include <iostream>
#include <iterator>
#include <exception>
#include <stdexcept>

#include "options_parser_normalize.h"
#include "options_normalize.h"

options_parser_normalize::options_parser_normalize() = default;



bool options_parser_normalize::parse( int argc, char ***argv, options *opts )
{
    options_normalize *opts_normalize = (options_normalize*) opts;

    namespace po = boost::program_options;
    po::variables_map vm;

    po::options_description desc( "PepSIRF: Peptide-based Serological Immune Response Framework score normalization module. \n"
                                );

    desc.add_options()
        ( "help,h", "Produce help message\n" );

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
