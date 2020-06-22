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

    po::options_description desc( "PepSIRF "
                                  + format_version_string()
                                  + ": Peptide-based Serological Immune Response "
                                  "Framework Paired (Duplicate) Enrichment module.\n",
                                  line_width
                                );
    desc.add_options()
        ( "help,h", "Produce help message and exit.\n"
          "The p_enrich module determines which peptides are enriched in samples that "
          "have been assayed in duplicate, as determined by user-specified thresholds. "
          "Thresholds are provided as comma-delimited pairs. In order for a peptide to "
          "be considered enriched, both replicates must meet or exceed the lower threshold "
          "and at least one replicate must meet or exceed the higher threshold, "
          "independent of order. Note that a peptide must meet each specified threshold "
          "(e.g., zscore, norm count and raw count) in order to be considered enriched.\n"
        )
        ( "samples,s", po::value( &opts_p_enrich->in_samples_fname )
          ->required(),
          "The name of the file containing sample pair information, denoting which "
          "samples, in the input matrices, are replicates. This file must be "
          "tab-delimited with one pair of samples per line.\n"
        )
        ( "zscores,z", po::value( &opts_p_enrich->in_zscore_fname )
          ->required(),
          "A tab-delimited matrix containing the zscores of each peptide in every "
          "sample of interest. This should be in the format output by the zscore module, "
          "with peptides on the rows and sample names on the columns.\n"
        )
        ( "zscore_constraint", po::value( &opts_p_enrich->zscore_params )
          ->required(),
          "Comma-separated zscores thresholds. These thresholds will be evaluated as "
          "specified in the module description. Example: '--zscore_constraint 3.5,4.0'.\n"
        )
        ( "norm_scores,n", po::value( &opts_p_enrich->in_norm_scores_fname )
          ->required(),
          "A tab-delimited matrix containing normalized counts for each peptide in each "
          "sample of interest.\n"
        )
        ( "norm_score_constraint", po::value( &opts_p_enrich->norm_scores_params )
          ->required(),
          "Comma-separated normalized count thresholds. These thresholds will be evaluated "
          "as specified in the module description. Example: '--norm_score_constraint 5.5,6.0'.\n"
        )
        ( "raw_scores,r", po::value( &opts_p_enrich->in_raw_scores_fname )
          ->default_value( "" ),
          "Optionally, a tab-delimited matrix containing raw counts can be included. This matrix "
          "must contain the raw counts for each peptide. If included, '--raw_score_constraint' "
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
          "The minimum total raw count a sample can have for all of its peptides in "
          "order for any of the peptides in that sample to be considered enriched. "
          "This provides a way to impose a minimum read count for a sample to be evaluated.\n"
        )
        ( "outfile_suffix,x",
          po::value( &opts_p_enrich->out_suffix )
          ->default_value( "" ),
          "Suffix to add to all output files. Together, the sample name and the suffix "
          "will form the name of the output file for each sample. For example, with a "
          "suffix of '_enriched.txt' and a sample name of ‘sample1’, the name of the output "
          "file for this sample would be ‘sample1_enriched.txt’. For example, "
          "'_enriched.txt' can be used. By default, no suffix is used.\n"
        )
        ( "join_on,j",
          po::value( &opts_p_enrich->out_fname_join )
          ->default_value( "~" ),
          "A character or string to use to join replicate sample names in order to create "
          "output file names. For a pair of samples, A and B, the resulting file will "
          "have the name 'A~B' if this flag is not given. Otherwise, the given value will "
          "be used in place of '~'.\n"
         )
        ( "output,o", po::value( &opts_p_enrich->out_dirname )
          ->default_value( "paired" ),
          "Directory name to which output files will be written. An output file will be "
          "generated for each sample with at least one enriched peptide. This directory "
          "will be created by the module.\n"
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
