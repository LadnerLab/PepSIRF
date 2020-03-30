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
         "Additional flags may be used to perform different analyses. "
         "For each input flag included, one output file will be written.\n"
        )
        ( "input,i", po::value( &opts_info->in_fname )->required(),
          "An input score matrix to gather information from.\n"
        )
        ( "get_samples,s", po::value( &opts_info->out_samples_fname )
          ->default_value( "" ),
          "Name of the file to write sample names to. Output will be "
          "in the form of a file with no header, one sample name per line.\n"
        )
        ( "get_probes,p", po::value( &opts_info->out_pep_names_fname ),
          "Name of the file to write probe names to. Output will be "
          "in the form of a file with no header, one probe name per line.\n"
        )
        ( "col_sums,c", po::value( &opts_info->out_col_sums_fname ),
          "Name of the file to write the sum of column scores to. "
          "Output will be a tab-delimited file with a header. The first "
          "entry in each column will be the name of the sample, and the second "
          "will be the sum of the scores each peptide had for the sample.\n"
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
