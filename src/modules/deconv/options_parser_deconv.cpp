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
        ( "help,h", "Produce help message\n" )
        ( "linked,l", po::value<std::string>( &opts_deconv->linked_fname ),
          "Name of file containing peptide to species linkages.\n"
        )
        ( "threshold,t", po::value<std::size_t>( &opts_deconv->threshold ),
          "Threshold number of peptides for a species to be considered.\n"
        )
        ( "output,o", po::value<std::string>( &opts_deconv->output_fname ),
          "Name of the file to write output to. Output will be in the form of "
          "a tab-delimited file with a header. Each entry will be of the form:\n"
          "species_id\\tcount\n"
        )
        ( "single_threaded", po::bool_switch( &opts_deconv->single_threaded )->default_value( false ),
          "By default this module uses two threads. Include this option with no arguments if you only want "
          " one thread to be used.\n"
        )
        ( "enriched,e", po::value<std::string>( &opts_deconv->enriched_fname ),
          "File containing the names of enriched peptides, one per line. "
          "Each file in this file should have a corresponding entry in the "
          "file provided by the --linked option.\n"
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

