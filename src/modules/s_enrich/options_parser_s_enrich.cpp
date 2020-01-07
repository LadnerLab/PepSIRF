#include <boost/program_options.hpp>
#include "options_parser_s_enrich.h"

bool options_parser_s_enrich
::parse( int argc, char ***argv, options *opts )
{
    options_s_enrich *opts_s_enrich = (options_s_enrich*) opts;
    namespace po = boost::program_options;
    po::variables_map vm;

    po::options_description desc( "PepSIRF: Peptide-based Serological Immune Response "
                                  "Framework species Single-Replicate Enrichment module"
                                );
    desc.add_options()
        ( "help,h", "Produce help message\n" )
        ( "zscores", po::value( &opts_s_enrich->in_zscore_fname )
          ->required(),
          ""
        )
        ( "min_zscore", po::value( &opts_s_enrich->min_zscore )
          ->required(),
          ""
        )
        ( "norm_scores", po::value( &opts_s_enrich->in_norm_score_fname )
          ->required(),
          ""
        )
        ( "min_norm_score", po::value( &opts_s_enrich->min_norm_score )
          ->required(),
          ""
        )
        ( "raw_counts", po::value( &opts_s_enrich->in_raw_count_fname )
          ->default_value( "" ),
          ""
        )
        ( "min_raw_count", po::value( &opts_s_enrich->min_raw_count )
          ->default_value( static_cast<double>( 0 ) )
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
