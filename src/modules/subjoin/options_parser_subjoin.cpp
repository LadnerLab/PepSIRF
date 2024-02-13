#include "options_parser_subjoin.h"
#include "file_io.h"
#include "logger.h"
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
        ( "multi_file,m", po::value<std::string>()->notifier( [&]( const std::string& provided_string )
                      {
                            if( !provided_string.empty() )
                                {
                                    std::ifstream input_f{ provided_string };
                                    if( input_f.fail() )
                                        {
                                            Log::error("Unable to open multi_file input.\n");
                                        }
                                    pepsirf_io::read_file(  input_f,
                                                            boost::is_any_of( "\t" ),
                                                            []( typename std::vector<std::string>::iterator a,
                                                                typename std::vector<std::string>::iterator end
                                                            )-> std::pair<std::string,std::string>
                                                            {
                                                                std::string f = *a;
                                                                ++a;
                                                                if( a != end )
                                                                    return std::make_pair( f, *a );
                                                                else
                                                                    return std::make_pair( f, "" );
                                                            },
                                                            std::back_inserter( opts_subjoin->multi_matrix_name_pairs )
                                                        );

                                }


                      }
                    ),
          "The name of a tab-delimited file containing score matrix "
          "and sample name list filename pairs, one per line. "
          "Each of these pairs must be a score matrix and a file "
          "containing the names of samples (or peptides, if specified) "
          "to keep from the input the score matrix. The score matrix should be of the format output by the "
          "demux module, with sample names on the columns and peptide names on the rows. "
          "The namelist must have one name per line, but can optionally have 2, if "
          "renaming samples in the subjoin output. Optionally, a name list can be "
          "omitted if all samples from the input matrix should be included in the "
          "output.\n"
        )
        ( "input,i", po::value( &matrix_name_list_pairs )->notifier(
                    [&]( const std::vector<std::string>& name_pairs )
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
                                              ->input_matrix_name_pairs.push_back( std::make_pair
                                                                             ( split_output[ 0 ],
                                                                               split_output[ 1 ]
                                                                               )
                                                                             );

                                      }
                                  else if( split_output.size() == 1 )
                                      {
                                          Log::warn(
                                            "No sample name list has been given for "
                                            + split_output[0]
                                            + ". All samples from this input matrix will be included.\n"
                                        );
                                          opts_subjoin->input_matrix_name_pairs.emplace_back( std::make_pair( split_output[0], "" ) );
                                      }
                                  else
                                      {
                                          throw po::validation_error(
                                          po::validation_error::invalid_option_value,
                                          "input",
                                          provided_value
                                          );

                                      }

                              }

                      }
                    ),
          "Comma-separated filenames (For example: score_matrix.tsv,sample_names.txt ). "
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
        ( "filter_peptide_names", po::bool_switch( &opts_subjoin->use_sample_names )->default_value( false )
          ->notifier( [&]( bool val ){ opts_subjoin->use_sample_names = !val; } ),
          "Flag to include if the name lists input to the input or multi_file options should be treated as "
          "peptide (i.e. row) names instead of sample (i.e. column) names. With the inclusion of this flag, the input files will "
          "be filtered on peptide names (rows) instead of sample names (column).\n"
        )
        (
         "exclude,e", po::bool_switch( &opts_subjoin->exclude_names)->default_value( false )
         ->notifier( [&]( bool val ){ opts_subjoin->exclude_names = !val; } ),
         "New data file will contain all of the input samples (or peptides) except the ones specified in the sample names.\n"
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
        ("logfile", po::value( &opts_subjoin->logfile )
         ->default_value( options_subjoin::set_default_log() ),
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
