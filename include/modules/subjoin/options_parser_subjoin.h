#ifndef OPTIONS_PARSER_SUBJOIN_HH_INCLUDED
#define OPTIONS_PARSER_SUBJOIN_HH_INCLUDED
#include <boost/program_options.hpp>

#include "options_parser.h"
#include "options_subjoin.h"

/**
 * Basic class for parsing options from the commandline. Makes use of 
 * boost::program_options. Note that this class is used for parsing options 
 * pertinent to the 'deconv' module.
 **/
class options_parser_subjoin : public options_parser
{
 public:

    /**
     * Default constructor.
     **/
    options_parser_subjoin();

    /**
     * Parse command-line arguments, store parsed options in the 
     * members of the opts item.
     * @param argc The number of arguments passed from the command-line.
     * @param argv Pointer to '**argv' obtained from the main function.
     * @param opts Reference to opts object, values will be populated by 
     *             the arguments from the commandline.
     * @throws error upon required argument not supplied.
     * @returns boolean true if arguments other than '-h, --help' have 
     *          been supplied, false otherwise.
     **/
    bool parse( int argc, char ***argv, options *opts );
};

#endif // OPTIONS_PARSER_SUBJOIN_HH_INCLUDED
