#include "options_parser_subjoin.h"
#include "file_io.h"
#include <boost/algorithm/string.hpp>
options_parser_subjoin::options_parser_subjoin() = default;


bool options_parser_subjoin::parse( int argc, char ***argv, options *opts )
{
    options_subjoin *opts_subjoin = (options_subjoin*) opts;
    namespace po = boost::program_options;
    po::variables_map vm;
    std::vector<std::string> matrix_name_list_pairs;

    po::options_description desc( "PepSIRF "
                                  + format_version_string()
                                  + ": Peptide-based Serological Immune "
                                  "Response Framework subjoin module",
                                  line_width
                                );

    desc.add_options()
        (
         "help,h", "Produce help message\n"
         "The subjoin module is used to manipulate matrix files. This module "
         "can create a subset of an existing matrix, can combine multiple matrices "
         "together or perform a combination of these two functions.\n"
        )
        ( "filter_scores,f", po::value( &matrix_name_list_pairs )->required()
          ->notifier( [&]( const std::vector<std::string>& name_pairs )
                      {
                          std::vector<std::string> split_output;

                          for( const auto& provided_value : name_pairs )
                              {
                                  boost::split( split_output,
                                                provided_value,
                                                boost::is_any_of( "," )
                                              );

                                  if( split_output.size() == 2 )
                                      {
                                          boost::trim( split_output[ 0 ] );
                                          boost::trim( split_output[ 1 ] );

                                          opts_subjoin
                                              ->matrix_name_pairs.push_back( std::make_pair
                                                                             ( split_output[ 0 ],
                                                                               split_output[ 1 ]
                                                                               )
                                                                             );

                                      }
                                  else if( split_output.size() == 1 )
                                      {
                                          std::ifstream input_f{ split_output[ 0 ] };
                                          if( input_f.fail() )
                                              {
                                                  throw std::runtime_error( "Unable to open input file for reading" );
                                              }

                                          pepsirf_io::read_file( input_f,
                                                                 boost::is_any_of( "\t" ),
                                                                 []( typename std::vector<std::string>::iterator a,
                                                                     typename std::vector<std::string>::iterator end
                                                                   )
                                                                 -> std::pair<std::string,std::string>
                                                                 { std::string f = *a;
                                                                     if( a != end )
                                                                         ++a;
                                                                   return std::make_pair( f, *a ); },
                                                                 std::back_inserter( opts_subjoin->matrix_name_pairs )
                                                               );
                                      }
                                  else
                                      {
                                          throw po::validation_error(
                                          po::validation_error::invalid_option_value,
                                          "filter_scores",
                                          provided_value
                                          );

                                      }

                              }

                      }
                    ),
          "Either comma-separated filenames (For example: score_matrix.tsv,sample_names.txt ), "
          "or the name of a tab-delimited file containing score_matrix "
          "and sample name list filename pairs, one per line. "
          "Each of these pairs must be a score matrix and a file "
          "containing the names of samples (or peptides, if specified) "
          "to keep in the score matrix. The score matrix should be of the format output by the "
          "demux module, with sample names on the columns and peptide names on the rows. "
          "The namelist must have one name per line, but can optionally have 2. If "
          "2 tab-delimited names are included on one line, the name in the first column "
          "should match the name in the input matrix file, while the name in the second "
          "column will be output. Therefore, this allows for the renaming of samples in "
          "the output. To use multiple name lists with multiple "
          "score matrices, include this argument multiple times. "
          "Optionally, a name list can be omitted if all samples from the input "
          "matrix should be included in the output.\n"
        )
        ( "filter_peptide_names", po::bool_switch( &opts_subjoin->use_sample_names  )->default_value( false )
          ->notifier( [&]( bool val ){ opts_subjoin->use_sample_names = !val; } ),
          "Flag to include if the name lists input to the filter_scores options should be treated as "
          "peptide (i.e. row) names instead of sample (i.e. column) names. With the inclusion of this flag, the input files will "
          "be filtered on peptide names (rows) instead of sample names (column).\n"
        )
        ( "duplicate_evaluation,d", po::value<std::string>()->default_value( "include" )
          ->notifier( [&]( const std::string& provided_value )
                      {
                          if( evaluation_strategy::is_valid( provided_value ) )
                              {
                                  opts_subjoin->duplicate_resolution_strategy =
                                      evaluation_strategy::from_string( provided_value );
                              }
                          else
                              {
                                  throw po::validation_error(
                                                             po::validation_error::invalid_option_value,
                                                             "duplicate_evaluation",
                                                             provided_value
                                                             );
                              }
                      }
                    ),
          "Defines what should be done when sample or peptide names are not unique across files being "
          "joined. Currently, three different duplicate evaluation strategies are available: \n"
          " - combine: Combine (with addition) the values associated with identical sample/peptide names "
          "from different files.\n\n"
          " - include: Include each duplicate, adding a suffix to the duplicate name detailing the "
          "file from which the sample came.\n\n"
          " - ignore: Ignore the possibility of duplicates. Behavior is undefined when duplicates are "
          " encountered in this mode Therefore, this mode is not recommended.\n\n"
        )
        (
         "output,o", po::value<std::string>( &opts_subjoin->out_matrix_fname )->default_value( "subjoin_output.tsv" ),
         "Name for the output score matrix file. The output will be in the form of the input, but with only the "
         "specified values (samplenames or peptides) found in the namelists. \n"
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
