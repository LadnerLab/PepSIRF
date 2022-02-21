#include <boost/program_options.hpp>
#include <iostream>
#include <iterator>
#include <exception>
#include <stdexcept>
#include <string>

#include "options_parser_qc.h"
#include "options_qc.h"

options_parser_qc::options_parser_qc() = default;

bool options_parser_qc::parse( int argc, char ***argv, options *opts )
{
    options_qc *opts_qc = (options_qc*) opts;

    namespace po = boost::program_options;
    po::variables_map vm;

    po::options_description desc(   "PepSIRF "
                                    + format_version_string()
                                    + ": Peptide-based Serological "
                                    "Immune Response Framework score qc module. \n", line_width
                                );

    desc.add_options()
        ( "help,h", "produce help message\n"
          ""
        )
        ( "fif,f", po::value<std::string>( &opts_qc->idx_fname ),
          "flex id filename"
          ""
        )
        ( "input_r1", po::value<std::string>( &opts_qc->input_r1_fname ),
          "first index filename input"
          ""
        )
        ( "input_r2", po::value<std::string>( &opts_qc->input_r2_fname ),
          "second index filename input"
          ""
        )
        ( "samplelist,s", po::value<std::string>( &opts_qc->samplelist_fname ),
          "samplelist filename"
          ""
        )
        ( "outfile,o", po::value<std::string>( &opts_qc->output_fname ),
          "output filename"
          ""
        )
        ( "sname", po::value<std::string>( &opts_qc->samplename )
                    ->default_value( "SampleName" ),
          "sample name"
          ""
        )
        ( "index1", po::value<std::string>()
                     ->default_value( "0,0,0" )
                     ->notifier( [&]( const std::string &vals ) {
                             opts_qc->set_info( &options_qc::index1_data,
                                                    vals
                                                 );
                                                                }
                               )
        )
        ( "index2", po::value<std::string>()
                     ->default_value( "0,0,0" )
                     ->notifier( [&]( const std::string &vals ) {
                             opts_qc->set_info( &options_qc::index2_data,
                                                    vals
                                                 );
                                                                }
                               )
        );


        po::store( po::command_line_parser( argc, *argv ).options( desc ).run(), vm);

        if( vm.count( "help" )
        || argc == 2 // argc == 2 when 'pepsirf qc' is called
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