#ifndef OPTIONS_PARSER_LINK2_HH_INCLUDED
#define OPTIONS_PARSER_LINK2_HH_INCLUDED
#include "options_parser.h"

class options_parser_link2 : public options_parser
{
 public:
    options_parser_link2();

    bool parse( int argc, char ***argv, options *opts );
    
};

#endif // OPTIONS_PARSER_LINK2_HH_INCLUDED
