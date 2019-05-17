#include "options_parser.h"

options_parser::options_parser() = default;

void options_parser::parse( int argc, char ***argv, options& opts )
{
    namespace po = boost::program_options;
    char **argv_loc = *argv;
    po::variables_map vm;

    po::options_description desc( "PepSIRF: Peptide-based Serological Immune Response Framework" );
    desc.add_options()
        ( "help", "Produce help message" )
        ( "input_r1", po::value<std::string>(), "Input forward reads fastq file to parse.")
        ( "input_r2", po::value<std::string>(), "Input reverse reads fastq file to parse.")
        ( "library,l", po::value<std::string>(), "Designed library containing amino acid peptides. "
                                                         "Library should be in fasta format and should contain "
                                                         "sequences that were used to design input_r1 and input_r2."
        )
        ( "num_threads,t", po::value<int>(), "Number of threads to use for analyses." );

    po::store( po::parse_command_line( argc, argv_loc, desc ), vm );

    if( vm.count( "help" ) )
        {
            std::cout << desc << std::endl;
        }
}
