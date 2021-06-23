#ifndef OPTIONS_PARSER_ENRICH_HH_INCLUDED
#define OPTIONS_PARSER_ENRICH_HH_INCLUDED
#include "options_parser.h"
#include "options_enrich.h"

class options_parser_enrich : public options_parser
{
    bool parse( int argc, char ***argv, options *opts );
};

#endif // OPTIONS_PARSER_ENRICH_HH_INCLUDED
