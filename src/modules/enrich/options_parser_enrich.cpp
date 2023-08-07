#include "options_parser_enrich.h"

#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>
#include <limits>
#include <stdexcept>
#include <iostream>
#include "stream_tools.h"
#include "file_io.h"
#include "predicate.h"

bool options_parser_enrich::parse( int argc, char ***argv, options *opts )
{
    options_enrich *opts_enrich = (options_enrich*) opts;

    namespace po = boost::program_options;
    po::variables_map vm;

    po::options_description desc( "PepSIRF "
                                  + format_version_string()
                                  + ": Peptide-based Serological Immune Response "
                                  "Framework Enrichment module.\n"
                                  "The enrich module determines which peptides are enriched in samples that "
                                  "have been assayed in n-replicate, as determined by user-specified thresholds. "
                                  "Thresholds can be provided either as single integers or comma-delimted integer pairs. "
                                  "In order for a peptide to be considered enriched, all replicates must meet or exceed the "
                                  "lower threshold and at least one replicate must meet or exceed the higher threshold, "
                                  "independent of order. Note that a peptide must meet each specified threshold "
                                  "(e.g., zscore, norm count, etc.) in order to be considered enriched. If any replicates "
                                  "fail the --raw_score_constraint, that sample will be omitted from the analysis.\n",
                                  line_width
                                );
    desc.add_options()
        ( "help,h", "Produce help message and exit.\n"
        )
        ("threshhold_file,t",
         po::value( &opts_enrich->threshold_fname )->required()
         ->notifier([&]( std::string input_filename )->void
             {
                std::ifstream input_f{ input_filename };
                if( input_f.fail() )
                {
                    Log::error("Unable to open threshold_file input.\n");
                }

                std::string line;
                std::vector<std::string> matrix_thresh_pairs;
                while( getline( input_f, line ) )
                {
                    boost::split( matrix_thresh_pairs, line, boost::is_any_of( "\t" ) );
                    if( matrix_thresh_pairs.size() == 2 )
                    {
                        opts_enrich->matrix_thresh_fname_pairs
                            .emplace_back(
                                std::make_pair(matrix_thresh_pairs[0], matrix_thresh_pairs[1]
                            ));
                    }
                    else
                    {
                        Log::error(
                            "The input file must contain a matrix filename"
                            " with threshold(s) where the filename and"
                            " thresholds are tab-delimited and thresholds are"
                            " comma-delimited.\n"
                        );
                    }
                }
             }),
          "The name of a tab-delimited file containing one tab-delimited matrix filename and threshold(s), one per line. "
          "If using more than one threshold for a given matrix, then separate by comma. A matrix file may contain any score of interest, as long as "
          "it includes scores for each peptide in every sample of interest, with peptides on the rows and sample names on the columns. "
          "The provided thresholds should be comma-separated if more than one is provided for a single matrix file.\n"
        )
        ( "samples,s", po::value( &opts_enrich->in_samples_fname ),
          "The name of the tab-delimited file containing sample information, denoting which "
          "samples, in the input matrices, are replicates. This file must be "
          "tab-delimited with each line containing a set of replicates.\n"
        )
        ( "raw_scores,r", po::value( &opts_enrich->in_raw_scores_fname )
          ->default_value( "" ),
          "Optionally, a tab-delimited matrix containing raw counts can be provided. This matrix "
          "must contain the raw counts for each Peptide for every sample of interest. If included, "
          "'--raw_score_constraint' must also be specified.\n"
        )
        ( "raw_score_constraint", po::value( &opts_enrich->raw_scores_params_str )
          ->default_value( "" )
          ->notifier( [&]( const std::string& params_str ) -> void
                      {
                          if( !predicate::biconditional( vm[ "raw_scores" ].defaulted(),
                                                         vm[ "raw_score_constraint" ].defaulted()
                                                       )

                              && params_str.empty() )
                              {
                                  throw std::runtime_error( "If either 'raw_scores' "
                                                            "or 'raw_score_constraint' options "
                                                            "are included, BOTH must be."
                                                            );
                              }
                      }
                    ),
          "The minimum total raw count across all peptides for a sample to be included in the analysis."
          "This provides a way to impose a minimum read count for a sample to be evaluated.\n"
        )
        ( "enrichment_failure_reason,f", po::value( &opts_enrich->out_enrichment_failure )
          ->default_value( "" ),
          "For each sample set that does not result in the generation of an enriched peptide file, "
          "a row of two tab-delimited columns is provided: the first column "
          "contains the replicate names (comma-delimited) and the second column provides the reason why "
          "the sample did not result in an enriched peptide file.\n"
          "This file is output to the same directory as the enriched peptide files. The 'Reason' column "
          "will contain one of the following: 'Raw read count threshold' or 'No enriched peptides'.\n"
        )
        ( "outfile_suffix,x",
          po::value( &opts_enrich->out_suffix )
          ->default_value( "" ),
          "Suffix to add to all output files. Together, the sample name(s) and the suffix "
          "will form the name of the output file for each sample. For example, with a "
          "suffix of '_enriched.txt' and a sample name of ‘sample1’, the name of the output "
          "file for this sample would be ‘sample1_enriched.txt’. "
          "By default, no suffix is used.\n"
        )
        ( "join_on,j",
          po::value( &opts_enrich->out_fname_join )
          ->default_value( "~" ),
          "A character or string to use to join replicate sample names in order to create "
          "output file names. For a pair of samples, 'A' and 'B', the resulting file will "
          "have the name 'A~B' if this flag is not given. Otherwise, the given value will "
          "be used in place of '~'.\n"
        )
        ( "output_filename_truncate", po::bool_switch( &opts_enrich->truncate_names )
          ->default_value( false ),
          "By default each filename in the output directory will include every replicate name "
          "joined by the 'join_on' value. Alternatively, if more than two replicates are "
          "being evaluated, then you may include this flag to stop the filenames from including "
          "more than 3 samplenames in the output. When this flag is used, the output names will "
          "be of the form 'A~B~C~1more', for example.\n"
        )
        ( "low_raw_reads,l", po::bool_switch( &opts_enrich->low_raw_reads )
          ->default_value( false ),
          "By default samples with any replicates below the raw read threshold will be dropped "
          "when this flag is included, replicates with reads above the threshold will be kept \n"
        )
        ( "output,o", po::value( &opts_enrich->out_dirname )
          ->default_value( "enriched" ),
          "Directory name to which output files will be written. An output file will be "
          "generated for each sample with at least one enriched peptide. This directory "
          "will be created by the module.\n"
        )
        ("logfile", po::value( &opts_enrich->logfile )
         ->default_value( options_enrich::set_default_log() ),
         "Designated file to which the module's processes are logged. By "
         "default, the logfile's name will include the module's name and the "
         "time the module started running.\n"
        )
        ;

    po::store( po::command_line_parser( argc, *argv ).options( desc ).run(), vm);

    if( vm.count( "help" )
        || argc == 2
        )
        {
            std::ostringstream info_str;
            info_str << desc << "\n";

            Log::info(info_str.str());

            return false;
        }
    else
        {
            po::notify( vm );
            return true;
        }
};
