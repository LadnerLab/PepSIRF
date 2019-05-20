#ifndef OPTIONS_PARSER_HH_INCLUDED
#define OPTIONS_PARSER_HH_INCLUDED
#include <boost/program_options.hpp>
#include <iostream>
#include <iterator>
#include <exception>
#include <stdexcept>

#include "options.h"

/*!
 * Basic class for parsing options from the commandline. Makes use of 
 * boost::program_options. 
 */
class options_parser
{
public:
    /**
     * Default constructor.
     **/
    options_parser(); 

    /**
     * Parse command-line arguments, store parsed options in the 
     * members of the opts item.
     * @param argc The number of arguments passed from the command-line.
     * @param argv Pointer to '**argv' obtained from the main function.
     * @param opts Reference to opts object, values will be populated by 
     *             the arguments from the commandline.
     * @throws error upon required argument not supplied.
     **/
    void parse( int argc, char ***argv, options& opts );
};





#endif //OPTIONS_PARSER_HH_INCLUDED
