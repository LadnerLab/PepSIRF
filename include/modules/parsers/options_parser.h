#ifndef OPTIONS_PARSER_HH_INCLUDED
#define OPTIONS_PARSER_HH_INCLUDED
#include <boost/program_options.hpp>
#include <iostream>
#include <iterator>
#include <exception>
#include <stdexcept>

#include "options.h"

/*!
 * Base class for parsing options from the commandline. Makes use of 
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
     * Default destructor.
     **/
    virtual ~options_parser();

    /**
     * Line width of the terminal.
     **/
    int line_width;

    /**
     * The current version of PepSIRF
     **/
    std::string version_no;

    /**
     * Parse command-line arguments, store parsed options in the 
     * members of the opts item.
     * @param argc The number of arguments passed from the command-line.
     * @param argv Pointer to '**argv' obtained from the main function.
     * @param opts Reference to opts object, values will be populated by 
     *             the arguments from the commandline.
     * @throws error upon required argument not supplied.
     * @returns boolean true if arguments other than '--help', '-h' are supplied, 
     *          false otherwise.
     **/
    virtual bool parse( int argc, char ***argv, options *opts );

    /**
     * Set the version number of PepSIRF.
     * @return a pointer to this
     **/
    options_parser *set_version( const std::string version_no );

    std::string format_version_string();
};





#endif //OPTIONS_PARSER_HH_INCLUDED
