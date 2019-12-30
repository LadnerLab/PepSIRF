#include "options_parser_zscore.h"
#include "options_zscore.h"

#include <boost/program_options.hpp>

options_parser_zscore::options_parser_zscore() = default;

bool options_parser_zscore::parse( int argc, char ***argv, options *opts )
{
    options_zscore *opts_zscore = (options_zscore*) opts;
    namespace po = boost::program_options;
    po::variables_map vm;
    std::vector<std::string> matrix_name_list_pairs;


    po::options_description desc( "PepSIRF: Peptide-based Serological Immune "
                                  "Response Framework zscore module"
                                );

    desc.add_options()
        (
         "help,h", "Produce help message\n"
        )
        ( "scores,s", po::value( &opts_zscore->in_fname )->required(),
          "Name of the file to use as input. Should be a score matrix in the "
          "format as output by the demux or subjoin modules.\n"
        )
        ( "bins,b", po::value( &opts_zscore->in_bins_fname )->required(),
          "Name of the file containing bins, one bin per line. "
          "Each bin contains a tab-separated list of peptide names.\n"
        )
        ( "output,o", po::value( &opts_zscore->out_fname )->default_value( "zscore_output.tsv" ),
          "Name of the file to write output to. In this file, "
          "each peptide will be written with its z-score within "
          "each sample.\n"
        )
        ( "num_threads,t", po::value( &opts_zscore->num_threads )
          ->default_value( opts_zscore->DEFAULT_NUM_THREADS ),
          "The number of threads to use for analyses.\n"
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
}
