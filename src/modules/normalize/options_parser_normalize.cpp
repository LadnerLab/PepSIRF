#include <boost/program_options.hpp>
#include <iostream>
#include <iterator>
#include <exception>
#include <stdexcept>
#include <string>

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
        ( "help,h", "Produce help message\n" )
        ( "peptide_scores", po::value<std::string>( &opts_normalize->peptide_scores_fname ),
          "Name of file containing peptide scores. This file should be tab-delimited "
          "with the first column being peptide names, and every next column should be \n"
          "the peptide's score within a given sample (the first item in the column). "
          "This is exactly the format output by the deconv module.\n"
        );

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
