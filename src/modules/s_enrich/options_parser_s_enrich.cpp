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

    po::options_description desc( "PepSIRF: Peptide-based Serological Immune Response "
                                  "Framework species Single-Replicate Enrichment module"
                                );
    desc.add_options()
        ( "help,h", "Produce help message and exit.\n"
          "This module determines which probes in a sample are "
          "enriched, as determined by various numerical thresholds that. "
          "Note that a probe must meet each specified numeric threshold in order "
          "to be considered enriched.\n"
        )
        ( "zscores", po::value( &opts_s_enrich->in_zscore_fname )
          ->required(),
          "A matrix containing the zscores of each probe in every sample. "
          "This should be in the format output by the zscore module, with "
          "probes on the rows and sample names on the columns.\n"
        )
        ( "min_zscore", po::value( &opts_s_enrich->min_zscore )
          ->required(),
          "The minimum zscore a probe must have in order to be considered "
          "enriched.\n"
        )
        ( "norm_scores", po::value( &opts_s_enrich->in_norm_score_fname )
          ->required(),
          "A matrix containing normalized scores for each probe in each sample.\n"
        )
        ( "min_norm_score", po::value( &opts_s_enrich->min_norm_score )
          ->required(),
          "The minimum normalized score a probe must have in a sample "
          "in order to be considered enriched.\n"
        )
        ( "raw_counts", po::value( &opts_s_enrich->in_raw_count_fname )
          ->default_value( "" ),
          "Optionally, a raw count matrix can be included. This matrix must "
          "contain the raw counts of each probe. If included, 'min_raw_count' "
          "must also be specified.\n"
        )
        ( "min_raw_count", po::value( &opts_s_enrich->min_raw_count )
          ->default_value( static_cast<double>( 0 ) )
          ->notifier( [&]( const double val ) -> void
                      {
                          if( !iff( vm[ "raw_counts" ].defaulted(),
                                    vm[ "min_raw_count" ].defaulted()
                                 )
                              // dirty hack to ensure the argument is used
                              && val < std::numeric_limits<double>::max()
                              )
                              {
                                  throw std::runtime_error( "If either 'raw_counts' "
                                                            "or 'min_raw_count' options "
                                                            "are included, BOTH must be."
                                                            );
                              }
                      }
                    ),
          "The minimum raw count a sample can have for all of its peptides in "
          "order for any of the probes in that sample to be considered enriched. "
          "The sum of each probe's raw count in a sample must be at least this "
          "value in order for the sample to be considered.\n"
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
