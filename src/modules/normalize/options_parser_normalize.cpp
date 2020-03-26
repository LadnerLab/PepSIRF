#include <boost/program_options.hpp>
#include <iostream>
#include <iterator>
#include <exception>
#include <stdexcept>
#include <string>

#include "options_parser_normalize.h"
#include "options_normalize.h"

options_parser_normalize::options_parser_normalize() = default;



bool options_parser_normalize::parse( int argc, char ***argv, options *opts )
{
    options_normalize *opts_normalize = (options_normalize*) opts;

    namespace po = boost::program_options;
    po::variables_map vm;

    po::options_description desc( "PepSIRF: Peptide-based Serological Immune Response Framework score normalization module. \n", line_width
                                );

    desc.add_options()
        ( "help,h", "Produce help message\n" )
        ( "peptide_scores,p", po::value<std::string>( &opts_normalize->peptide_scores_fname ),
          "Name of file containing peptide scores. This file should be tab-delimited "
          "with the first column being peptide names, and every next column should be \n"
          "the peptide's score within a given sample (the first item in the column). "
          "This is exactly the format output by the deconv module.\n"
        )
        ( "col_sum,c", po::bool_switch( &opts_normalize->col_sum_norm )->default_value( true ),
          "Normalize the counts using a column-sum method. Output is the number of reads a peptide "
          "per reads million mapped. Note that if size_factors is also included, the value of this flag "
          "will be ignored and the size_factors method is used. By default, col_sum normalization is used.\n"
        )
        ( "size_factors,s", po::bool_switch( &opts_normalize->size_factors_norm )->default_value( false ),
          "Normalize the counts using the size factors method (Anders and Huber 2010). Note that if this "
          "flag is included, the value of col_sum will be ignored.\n"
        )
        ( "precision,p", po::value( &opts_normalize->precision_digits )
          ->default_value( static_cast<std::size_t>( 2 ) ),
          "Output score precision. The scores written to the output will be output to this "
          "many decimal places. For example, a value of 2 will result in values output in the form "
          "'23.45', and a value of 3 will result in output of the form '23.449.\n"
        )
        ( "output,o", po::value<std::string>( &opts_normalize->output_fname )
          ->default_value( "norm_output.tsv" ),
          "The name of the file to write output to. The output is formatted in the same "
          "way the input 'peptide_scores' are formatted, i.e. a score matrix with samples "
          "on the columns and scores for a certain peptide on the rows. The score for each peptide "
          "in the output has been normalized in the manner specified.\n"
        );

    po::store( po::command_line_parser( argc, *argv ).options( desc ).run(), vm);

    if( vm.count( "help" )
        || argc == 2 // argc == 2 when 'pepsirf norm' is called
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
