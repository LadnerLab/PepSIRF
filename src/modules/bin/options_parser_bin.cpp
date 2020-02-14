#include "options_parser_bin.h"
#include "options_bin.h"
#include <boost/program_options.hpp>


options_parser_bin::options_parser_bin() = default;


bool options_parser_bin::parse( int argc, char ***argv, options *opts )
{
    options_bin *opts_bin = (options_bin*) opts;
    namespace po = boost::program_options;
    po::variables_map vm;

    po::options_description desc( "PepSIRF: Peptide-based Serological Immune "
                                  "Response Framework bin module.\n",
                                  line_width
                                );

    desc.add_options()
        (
         "help,h", "Produce help message and exit.\n"
         "Typically, the subjoin module will be used to "
         "specify samples that are negative controls. \n"
         "This module is usually used to create bins from these, "
         "negative controls, allowing the z-scores for a peptide \n"
         "to be calculated from the scores of the peptides it shares a bin with.\n"

        )
        (
         "scores,s", po::value( &opts_bin->input_scores_fname )->required(),
         "The input file to parse scores from. Bins will be created from the peptides "
         "in these bins and the scores the peptides have. Peptides with similar scores will be "
         "binned together.\n"
        )
        (
         "bin_size,b", po::value( &opts_bin->min_bin_size )
         ->default_value( opts_bin->DEFAULT_MIN_BIN_SIZE ),
         "The minimum number of peptides that can be placed in a bin. "
         "If a bin would be created with fewer than bin_size peptides, "
         "it will be combined with the next bin until at least bin_size "
         "peptides are found.\n"
        )
        (
         "round_to,r", po::value( &opts_bin->rounding_factor )
         ->default_value( opts_bin->DEFAULT_ROUNDING_FACTOR ),
         "The 'rounding factor' for the scores parsed from the score matrix. "
         "Scores found in the matrix will be rounded to the nearest 1/10^x "
         "for a rounding factor x. For example, a rounding factor of 0 "
         "will result in rounding to the nearest integer, while a rounding "
         "factor of 1 will result in rounding to the nearest tenth.\n"
        )
        (
          "output,o", po::value( &opts_bin->output_bins_fname )->default_value( "bin_output.tsv" ),
          "Name of the file to write the bins to, one per line. Each bin "
          "will be a tab-delimited list of the names of the peptides in the bin.\n"
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
