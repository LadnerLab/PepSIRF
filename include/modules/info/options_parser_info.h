#ifndef OPTIONS_PARSER_INFO_HH_INCLUDED
#define OPTIONS_PARSER_INFO_HH_INCLUDED
#include <boost/program_options.hpp>

#include "options_info.h"
#include "options_parser.h"

class options_parser_info : public options_parser
{

public:
    options_parser_info() = default;

    bool parse( int argc, char ***argv, options *opts );

};

#endif // OPTIONS_PARSER_INFO_HH_INCLUDED
