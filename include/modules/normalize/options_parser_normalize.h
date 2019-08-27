#ifndef OPTIONS_PARSER_NORMALIZE_HH_INCLUDED
#define OPTIONS_PARSER_NORMALIZE_HH_INCLUDED
#include "options_parser.h"
#include "options.h"

class options_parser_normalize : public options_parser
{
 public:
    options_parser_normalize();

    bool parse( int argc, char ***argv, options *opts );


};

#endif // OPTIONS_PARSER_NORMALIZE_HH_INCLUDED
