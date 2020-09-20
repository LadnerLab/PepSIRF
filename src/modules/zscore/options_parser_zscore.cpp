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


    po::options_description desc( "PepSIRF "
                                  + format_version_string()
                                  + ": Peptide-based Serological Immune "
                                  "Response Framework zscore module",
                                  line_width
                                );

    desc.add_options()
        (
         "help,h", "Produce help message\n"
         "The zscore module is used to calculate Z scores for each peptide in each sample. "
         "These Z scores represent the number of standard deviations away from the mean, "
         "with the mean and standard deviation both calculated separately for each bin of peptides.\n"
        )
        ( "scores,s", po::value( &opts_zscore->in_fname )->required(),
          "Name of the file to use as input. Should be a score matrix in the "
          "format as output by the demux and subjoin modules. "
          "Raw or normalized read counts can be used.\n"
        )
        ( "bins,b", po::value( &opts_zscore->in_bins_fname )->required(),
          "Name of the file containing bins, one bin per line, as output by the bin module. "
          "Each bin contains a tab-delimited list of peptide names.\n"
        )
        ( "trim,t", po::value( &opts_zscore->trim_percent )
          ->default_value( 2.5 )
          ->notifier( []( const double& val )
                      {
                          if( !( val >= 0.0 && val <= 100.00 ) )
                              {

                                  throw po::validation_error(
                                    po::validation_error::invalid_option_value,
                                    "filter_scores",
                                    std::to_string( val )
                                                             );
                              }

                      }
                      ),
          "Percentile of lowest and highest counts within a bin to ignore "
          "when calculating the mean and standard deviation. This value must be "
          "in the range [0.00,100.0].\n"
        )
        ("hdi,d", po::value( &opts_zscore->hpd_percent )->default_value( 0.0 )
          ->notifier( []( const double& val )
                      {
                          if( !( val >= 0.0 && val <= 100.00 ) )
                              {

                                  throw po::validation_error(
                                    po::validation_error::invalid_option_value,
                                    "filter_scores",
                                    std::to_string( val )
                                                             );
                              }

                      }
                      ),
          "Alternative approach for discarding outliers prior to calculating mean and stdev. "
          "If provided, this argument will override --trim, which trims evenly from both sides "
          "of the distribution. For --hdi, the user should provide the high density interval to "
          "be used for calculation of mean and stdev. For example, \"--hdi 0.95\" would instruct "
          "the program to utilize the 95% highest density interval (from each bin) for these calculations.\n"
        )
        ( "output,o", po::value( &opts_zscore->out_fname )->default_value( "zscore_output.tsv" ),
          "Name for the output Z scores file. This file will be a tab-delimited matrix file with "
          "the same dimensions as the input score file. Each peptide will be written with its "
          "z-score within each sample.\n"
        )
        (
         "nan_report,n", po::value( &opts_zscore->nan_report_fname )->default_value( "" ),
         "Name of the file to write out information regarding peptides that are given a zscore of 'nan'. "
         "This occurs when the mean score of a bin and the score of the focal peptide are both zero. "
         "This will be a tab-delimited file, with three columns per line. The first column will "
         "contain the name of the peptide, the second will be the name of the sample, and the third "
         "will be the bin number of the probe. This bin number corresponds to the line number in the "
         "bins file, within which the probe was found.\n"
        )
        ( "num_threads", po::value( &opts_zscore->num_threads )
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
