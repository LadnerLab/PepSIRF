#ifndef OPTIONS_PARSER_NORMALIZE_HH_INCLUDED
#define OPTIONS_PARSER_NORMALIZE_HH_INCLUDED
#include "options_parser.h"
#include "options.h"

/**
 * Options parser for the normalization module.
 **/
class options_parser_normalize : public options_parser
{
 public:

    /**
     * Default constructor.
     **/
    options_parser_normalize();

    /**
     * Parse the options provided by the command-line, which are 
     * stored in argc/argv. Values are stored in opts. 
     * @param argc The number of args.
     * @param argv Pointer to argv.
     * @param opts Pointer to options object 
     * @returns boolean true if normalization should happen, false otherwise.
     **/
    bool parse( int argc, char ***argv, options *opts );


};

#endif // OPTIONS_PARSER_NORMALIZE_HH_INCLUDED
