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
        ( "index_fname,i", "flex id filename"
          ""
        )
        ( "input_r1", "first index filename input"
          ""
        )
        ( "input_r2", "second index filename input"
          ""
        )
        ( "samplelist,s", "samplelist filename"
          ""
        )
        ( "output,o", "output filename"
          ""
        );

}