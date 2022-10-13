#include "options_parser_info.h"
#include <boost/program_options.hpp>

options_parser_info::options_parser_info() = default;
options_parser_info::~options_parser_info() = default;

bool options_parser_info
::parse( int argc, char ***argv, options *opts )
{
    options_info *opts_info = (options_info*) opts;
    namespace po = boost::program_options;
    po::variables_map vm;

    po::options_description desc( "PepSIRF "
                                  + format_version_string()
                                  + ": Peptide-based Serological Immune "
                                  "Response Framework info module",
                                  line_width
                                );

    desc.add_options()
        (
         "help,h", "Produce help message and exit.\n"
         "This module is used to gather information about a score matrix. "
         "By default, the number of samples and peptides in the matrix will be output. "
         "Additional flags may be provided to extract different types of information. "
         "Each of these flags should be accompanied by an output file name, "
         "to which the information be written.\n"
        )
        ( "input,i", po::value( &opts_info->in_fname )->required(),
          "An input score matrix to gather information from.\n"
        )
        ( "get_samples,s", po::value( &opts_info->out_samples_fname )
          ->default_value( "" ),
          "Name of a file to which sample names (i.e., column headers) should be written. Output will be "
          "in the form of a file with no header, one sample name per line.\n"
        )
        ( "get_probes,p", po::value( &opts_info->out_pep_names_fname ),
          "Name of a file to which peptide/probe names (i.e., row names) should be written. "
          "Output will be in the form of a file with no header, one peptide/probe name per line.\n"
        )
        ( "col_sums,c", po::value( &opts_info->out_col_sums_fname ),
          "Name of a file to which the sum of column scores should be written. "
          "Output will be a tab-delimited file with a header. The first "
          "entry in each column will be the name of the sample, and the second "
          "will be the sum of the peptide/probe scores for the sample.\n"
        )
        ("rep_names,n", po::value(&opts_info->in_replicates_fname),
            "An input file that the sequence names of replicates can be found. "
            "This file is required to run -a, --get_avgs. \n"
            )
        ( "get_avgs,a", po::value( &opts_info->out_avgs_fname ),
          "Name of a file to which the average of different replicate values should be written. "
          "Output will be a tab-delimted file with sample names as the column headers and "
          "peptide names as row names. Each entry consists of the average of the replicate "
          "values for the given sample and peptide. \n"
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

    return true;
}
