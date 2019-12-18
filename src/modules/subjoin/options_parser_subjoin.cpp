#include "options_parser_subjoin.h"
#include <boost/algorithm/string.hpp>
options_parser_subjoin::options_parser_subjoin() = default;


bool options_parser_subjoin::parse( int argc, char ***argv, options *opts )
{
    options_subjoin *opts_subjoin = (options_subjoin*) opts;
    namespace po = boost::program_options;
    po::variables_map vm;
    std::vector<std::string> matrix_name_list_pairs;


    po::options_description desc( "PepSIRF: Peptide-based Serological Immune "
                                  "Response Framework subjoin module"
                                );

    desc.add_options()
        (
         "help,h", "Produce help message\n"
        )
        ( "filter_scores", po::value( &matrix_name_list_pairs )->required()
          ->notifier( [&]( const std::vector<std::string>& name_pairs )
                      {
                          std::vector<std::string> split_output;

                          for( const auto& provided_value : name_pairs )
                              {
                                  boost::split( split_output,
                                                provided_value,
                                                boost::is_any_of( "," )
                                              );

                                  if( split_output.size() != 2 )
                                      {
                                          std::cout << split_output.size() << std::endl;
                                          throw po::validation_error(
                                          po::validation_error::invalid_option_value,
                                          "filter_scores",
                                          provided_value
                                                                     );                                         
                                      }

                                  boost::trim( split_output[ 0 ] );
                                  boost::trim( split_output[ 1 ] );

                                  opts_subjoin
                                  ->matrix_name_pairs.push_back( std::make_pair
                                                                 ( split_output[ 0 ],
                                                                   split_output[ 1 ]
                                                                 )
                                                               );
                              }
                          
                      }
                    ),
          "Comma-separated filenames (For example: score_matrix.tsv,peptide_names.txt ) "
          "for a score matrix and a file containing the names of peptides "
          "to keep in the score matrix. The score matrix should be of the format output by the "
          "demux module, with sample names on the columns and peptide names on the rows. "
          "The peptide namelist must have one name per line. To use multiple name lists with multiple "
          "score matrices, include this argument multiple times.\n"
        )
        (
         "output", po::value<std::string>( &opts_subjoin->out_matrix_fname )->default_value( "subjoin_output.tsv" ),
         "The name of the file to write output scores to. The output will be in the form of the input, but with only the "
         "peptides found in the namelists. \n"
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
