#ifndef OPTIONS_PARSER_LINK_HH_INCLUDED
#define OPTIONS_PARSER_LINK_HH_INCLUDED
#include "options_parser.h"

class options_parser_link : public options_parser
{
 public:
    options_parser_link();

    bool parse( int argc, char ***argv, options *opts );
    
};

#endif // OPTIONS_PARSER_LINK_HH_INCLUDED
