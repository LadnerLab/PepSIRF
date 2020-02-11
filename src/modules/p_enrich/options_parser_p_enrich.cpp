#include "options_parser_p_enrich.h"

#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>
#include <limits>
#include <stdexcept>
#include <iostream>
#include "stream_tools.h"

#include "predicate.h"

bool options_parser_p_enrich::parse( int argc, char ***argv, options *opts )
{
    options_p_enrich *opts_p_enrich = (options_p_enrich*) opts;

    namespace po = boost::program_options;
    po::variables_map vm;

    po::options_description desc( "PepSIRF: Peptide-based Serological Immune Response "
                                  "Framework species Paired (Duplicate) Enrichment module. "
                                  "This module determines which probes in samples done in duplicate are "
                                  "enriched, as determined by various numerical thresholds. "
                                  "For this module, thresholds are given as a comma-delimited pair. In order for the threshold "
                                  "to be met, each individual threshold in the pair must be met by at least "
                                  "one of the probes under consideration, independent of order. "
                                  "Note that a probe must meet each specified numeric threshold in order "
                                  "to be considered enriched.\n"
                                );
    desc.add_options()
        ( "help,h", "Produce help message and exit.\n" )
        ( "samples,s", po::value( &opts_p_enrich->in_samples_fname )
          ->required(),
          "The name of the file containing sample pairs, denoting which "
          "samples are in duplicate. This file must be tab-delimited with "
          "one pair or samples per line.\n"
        )
        ( "zscores,z", po::value( &opts_p_enrich->in_zscore_fname )
          ->required(),
          "A matrix containing the zscores of each probe in every sample. "
          "This should be in the format output by the zscore module, with "
          "probes on the rows and sample names on the columns.\n"
        )
        ( "zscore_constraint", po::value( &opts_p_enrich->zscore_params )
          ->required(),
          "Comma-separated values zscores should be constrained by. "
          "These scores will be evaluated as specified in the module description. "
          "Example: '--zscore_constraint 3.5,4.0'. \n"
        )
        ( "norm_scores,n", po::value( &opts_p_enrich->in_norm_scores_fname )
          ->required(),
          "A matrix containing normalized scores for each probe in each sample.\n"
        )
        ( "norm_score_constraint", po::value( &opts_p_enrich->norm_scores_params )
          ->required(),
          "Comma-separated values normalized scores should be constrained by. "
          "These scores will be evaluated as specified in the module description. "
          "Example: '--norm_score_constraint 5.5,6.0'. \n"
        )
        ( "raw_scores,r", po::value( &opts_p_enrich->in_raw_scores_fname )
          ->default_value( "" ),
          "Optionally, a raw count matrix can be included. This matrix must "
          "contain the raw counts of each probe. If included, 'min_raw_count' "
          "must also be specified.\n"
        )
        ( "raw_score_constraint", po::value( &opts_p_enrich->raw_scores_params )
          ->default_value( std::pair<double,double>{ 0.0, 0.0 } )
          ->notifier( [&]( const std::pair<double,double>& val ) -> void
                      {
                          if( !predicate::biconditional( vm[ "raw_scores" ].defaulted(),
                                                         vm[ "raw_score_constraint" ].defaulted()
                                                       )

                              // dirty hack to ensure the argument is used
                              && val.first < std::numeric_limits<double>::max()
                              )
                              {
                                  throw std::runtime_error( "If either 'raw_scores' "
                                                            "or 'raw_score_constraint' options "
                                                            "are included, BOTH must be."
                                                            );
                              }
                      }
                    ),
          "The minimum raw count a sample can have for all of its peptides in "
          "order for any of the probes in that sample to be considered enriched. "
          "The sum of each probe's raw count in each sample must be at least either of these "
          "values in order for the sample to be considered.\n"
        )
        ( "outfile_suffix,s",
          po::value( &opts_p_enrich->out_suffix )
          ->default_value( "" ),
          "Suffix to add to the names of the samples "
          "written to output. For example, '_enriched.txt' can be used. "
          "By default, no suffix is used.\n"
        ) 
        ( "join_on,j",
          po::value( &opts_p_enrich->out_fname_join )
          ->default_value( "~" ),
          "A character or string to join output sample names on. "
          "For a sample pair of samples A and B, the resulting file will "
          "have the name 'A~B' if this flag is not given. Otherwise, "
          "the given value will be used.\n"
         )
        ( "output,o", po::value( &opts_p_enrich->out_dirname )
          ->default_value( "paired" ),
          "Name of the directory to write output files to. "
          "Each sample with at least one enriched peptide will "
          "receive a file in the output directory.\n"
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
};
