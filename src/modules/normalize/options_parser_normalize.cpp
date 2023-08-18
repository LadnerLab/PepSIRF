#include <boost/program_options.hpp>
#include <iostream>
#include <iterator>
#include <exception>
#include <stdexcept>
#include <string>

#include "logger.h"
#include "options_parser_normalize.h"
#include "options_normalize.h"

options_parser_normalize::options_parser_normalize() = default;

bool options_parser_normalize::parse( int argc, char ***argv, options *opts )
{
    options_normalize *opts_normalize = (options_normalize*) opts;

    namespace po = boost::program_options;
    po::variables_map vm;

    po::options_description desc( "PepSIRF "
                                  + format_version_string()
                                  + ": Peptide-based Serological "
                                  "Immune Response Framework score normalization module. \n", line_width
                                );

    desc.add_options()
        ( "help,h", "Produce help message\n"
          "The norm module is used to normalize raw count data to allow for meaningful "
          "comparison among samples.\n"
        )
        ( "peptide_scores,p", po::value<std::string>( &opts_normalize->peptide_scores_fname ),
          "Name of tab-delimited matrix file containing peptide scores. "
          "This file should be in the same format as the output from the demux module.\n"
        )
        ( "normalize_approach,a", po::value<std::string>( &opts_normalize->approach )->default_value( "col_sum" ),
          "Options: \'col_sum\', \'size_factors\', \'diff\', \'ratio\', \'diff_ratio\', By default, col_sum normalization is used.\n"
          "\'col_sum\': Normalize the scores using a column-sum method. Output per peptide is the score per million for the sample "
          "(i.e., summed across all peptides).\n"
          "\'size_factors\': Normalize the scores using the size factors method (Anders and Huber 2010).\n"
          "\'diff\': Normalize the scores using the difference method. For each peptide and sample, the "
          "difference between the score and the respective peptide\'s mean score in the negative controls is determined.\n"
          "\'ratio\': Normalize the scores using the ratio method. For each peptide and sample, the ratio of score to the respective "
          "peptide\'s mean score in the negative controls is determined.\n"
          "\'diff_ratio\': Normalize the scores using the difference-ratio method. For each peptide and sample, the difference between "
          "the score and the respective peptide\'s mean score in the negative controls is first determined. This difference is then "
          "divided by the respective peptide\'s mean score in the negative controls.\n"
        )
        ( "negative_control", po::value<std::string>( &opts_normalize->neg_control )->default_value( "" ),
          "Optional data matrix for sb samples.\n"
        )
        ( "negative_id,s", po::value<std::string>( &opts_normalize->neg_id )->default_value( "" ),
          "Optional approach for identifying negative controls. Provide a unique string at the start of all negative control samples.\n"
        )
        ( "negative_names,n", po::value<std::string>( &opts_normalize->neg_names )->default_value( "" ),
          "Optional approach for identifying negative controls. Comma-separated list of negative control sample names.\n"
        )
        ( "precision", po::value( &opts_normalize->precision_digits )
          ->default_value( static_cast<std::size_t>( 2 ) ),
          "Output score precision. The scores written to the output will be output to this "
          "many decimal places. For example, a value of 2 will result in values output in the form "
          "'23.45', and a value of 3 will result in output of the form '23.449.\n"
        )
        ( "output,o", po::value<std::string>( &opts_normalize->output_fname )
          ->default_value( "norm_output.tsv" ),
          "Name for the output file. The output is formatted in the same "
          "way the input file provided with 'peptide_scores' (i.e., a score matrix with samples "
          "on the columns and scores for a certain peptide on the rows). The score for each peptide "
          "in the output has been normalized in the manner specified.\n"
        )
        ( "logfile", po::value( &opts_normalize->logfile )
          ->default_value( options_normalize::set_default_log() ),
          "Designated file to which the module's processes are logged. By "
          "default, the logfile's name will include the module's name and the "
          "time the module started running.\n"
        )
        ;

    po::store( po::command_line_parser( argc, *argv ).options( desc ).run(), vm);

    if( vm.count( "help" )
        || argc == 2 // argc == 2 when 'pepsirf norm' is called
      )
        {
            std::ostringstream info_str;
            info_str << desc << std::endl;

            Log::info(info_str.str());

            return false;
        }
    else
        {
            po::notify( vm );
            return true;
        }
}
