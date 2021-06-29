#include <boost/program_options.hpp>
#include <limits>
#include <stdexcept>
#include "options_parser_s_enrich.h"

bool options_parser_s_enrich
::parse( int argc, char ***argv, options *opts )
{
    // biconditional p <==> q
    auto iff = [] ( const bool p, const bool q )
        {
            return ( !p | q ) && ( !q | p );
        };

    options_s_enrich *opts_s_enrich = (options_s_enrich*) opts;
    namespace po = boost::program_options;
    po::variables_map vm;

    po::options_description desc( "PepSIRF "
                                  + format_version_string()
                                  + ": Peptide-based Serological Immune Response "
                                  "Framework species Single-Replicate Enrichment module",
                                  line_width
                                );
    desc.add_options()
        ( "help,h", "Produce help message and exit.\n"
          "The s_enrich module determines which peptides are enriched in each sample, "
          "as determined by user-specified thresholds. "
          "Note that a peptide must meet each specified threshold (e.g., zscore, "
          "norm count and raw count) in order to be considered enriched. This module "
          "will batch process all samples contained within a matrix and will generate "
          "one output file per sample.\n"
        )
        ( "zscores,z", po::value( &opts_s_enrich->in_zscore_fname )
          ->required(),
          "A tab-delimited matrix containing the zscores of each peptide in every "
          "sample of interest. This should be in the format output by the zscore module, "
          "with peptides on the rows and sample names on the columns.\n"
        )
        ( "min_zscore", po::value( &opts_s_enrich->min_zscore )
          ->required(),
          "The minimum zscore a peptide must have in order to be considered "
          "enriched.\n"
        )
        ( "norm_scores,n", po::value( &opts_s_enrich->in_norm_score_fname )
          ->required(),
          "A tab-delimited matrix containing normalized counts for each peptide in each "
          "sample of interest.\n"
        )
        ( "min_norm_score", po::value( &opts_s_enrich->min_norm_score )
          ->required(),
          "The minimum normalized count a peptide must have in a sample "
          "in order to be considered enriched.\n"
        )
        ( "raw_scores,r", po::value( &opts_s_enrich->in_raw_count_fname )
          ->default_value( "" ),
          "Optionally, a tab-delimited matrix containing raw counts can be included. "
          "This matrix must contain the raw counts for each peptide. If included, "
          "'min_raw_score' must also be specified.\n"
        )
        ( "min_raw_score", po::value( &opts_s_enrich->min_raw_count )
          ->default_value( static_cast<double>( 0 ) )
          ->notifier( [&]( const double val ) -> void
                      {
                          if( !iff( vm[ "raw_scores" ].defaulted(),
                                    vm[ "min_raw_score" ].defaulted()
                                 )
                              // dirty hack to ensure the argument is used
                              && val < std::numeric_limits<double>::max()
                              )
                              {
                                  throw std::runtime_error( "If either 'raw_scores' "
                                                            "or 'min_raw_score' options "
                                                            "are included, BOTH must be."
                                                            );
                              }
                      }
                    ),
          "The minimum total raw count a sample can have for all of its peptides in "
          "order for any of the peptides in that sample to be considered enriched. "
          "This provides a way to impose a minimum read count for a sample to be "
          "evaluated.\n"
        )
        ( "outfile_suffix,s",
          po::value( &opts_s_enrich->out_suffix )
          ->default_value( "" ),
          "Suffix to add to all output files. Together, the sample name and the suffix "
          "will form the name of the output file for each sample. For example, with a suffix "
          "of '_enriched.txt' and a sample name of 'sample1', the name of the output "
          "file for this sample would be 'sample1_enriched.txt. "
          "By default, no suffix is used.\n"
        )
        ( "output,o", po::value( &opts_s_enrich->out_dirname )
          ->default_value( "single" ),
          "Directory name to which output files will be written. An output file will be "
          "generated for each sample with at least one enriched peptide. "
          "This directory will be created by the module.\n"
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
