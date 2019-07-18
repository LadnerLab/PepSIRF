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
        ( "threshold,t", po::value<float>( &opts_deconv->threshold ),
          "Threshold number of peptides for a species to be considered.\n"
        )
        ( "output,o", po::value<std::string>( &opts_deconv->output_fname )->default_value( "deconv_output.tsv" ),
          "Name of the file to write output to. Output will be in the form of "
          "a tab-delimited file with a header. Each entry will be of the form:\n"
          "species_id\\tcount\n"
        )
        ( "single_threaded", po::bool_switch( &opts_deconv->single_threaded )->default_value( false ),
          "By default this module uses two threads. Include this option with no arguments if you only want "
          " one thread to be used.\n"
        )
        ( "fractional_scoring", po::bool_switch( &opts_deconv->fractional_scoring )->default_value( false ),
          "Use fractional instead of integer scoring. For integer scoring the score of each species is "
          "defined by the number of peptides that share a 7mer with that species. For fractional scoring "
          "the score of each species is defined by 1/n for each peptide, where n is the number of species "
          "a peptide shares a 7mer with. In this method of scoring "
          "peptides with fewer species are worth more.\n" 
        )
        ( "summation_scoring", po::bool_switch( &opts_deconv->summation_scoring )->default_value( false ),
          ""
        )
        ( "enriched,e", po::value<std::string>( &opts_deconv->enriched_fname ),
          "File containing the names of enriched peptides, one per line. "
          "Each file in this file should have a corresponding entry in the "
          "file provided by the --linked option.\n"
        )
        ( "create_linkage", po::bool_switch( &opts_deconv->create_linkage )->default_value( false ),
          "Boolean switch to create the linkage file that is used as input for "
          "the species deconvolution process. If this option is included "
          "then 'protein_file' and 'peptide_file' must also be included.\n"
        )
        ( "protein_file", po::value<std::string>( &opts_deconv->prot_file_fname ),
          "Name of fasta file containing protein sequences from which a design was "
          "created.\n"
        )
        ( "peptide_file", po::value<std::string>( &opts_deconv->peptide_file_fname ),
          "Name of fasta file containing aa peptides that have been designed as part "
          "of a library.\n"
        )
        ( "k_size,k", po::value<std::size_t>( &opts_deconv->k ), "Kmer size to use.\n" )
        ( "id_name_map", po::value<std::string>( &opts_deconv->id_name_map_fname )->default_value( "" ),
          "Map with tab-delimited lines, where the first item in each line is a species and the second is the name "
          "of the species.\n"
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

