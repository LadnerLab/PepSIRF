#ifndef OPTIONS_PARSER_FACTORY_HH_INCLUDED 
#define OPTIONS_PARSER_FACTORY_HH_INCLUDED 

#include <boost/program_options.hpp>
#include <algorithm>
#include <string>
#include <stdexcept>

#include "options_parser.h"
#include "options_parser_demux.h"

/**
 * Factory for creating parsers. Given commandline arguments passed
 * creates a factory that will create a commandline arguments object based 
 * upon the module passed in.
 **/
class options_parser_factory
{
 public:
    /**
     * Default constructor.
     **/
    options_parser_factory();

    /**
     * Create a parser for the module specified by the commandline.
     * Parses the name of the requested module from argv, creates and returns 
     * parser of the appropriate type.
     * @arg argc The number of arguments passed from the command-line.
     * @arg argv Pointer to argv, which is passed to the 'main' function.
     * @throws Exception if no appropriate module type is passed in from the commandline. 
     **/
    options_parser *create( int argc, char ***argv );
};

#endif // OPTIONS_FACTORY_HH_INCLUDED 
