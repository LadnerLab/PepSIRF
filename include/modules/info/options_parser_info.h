#ifndef OPTIONS_PARSER_INFO_HH_INCLUDED
#define OPTIONS_PARSER_INFO_HH_INCLUDED

#include "options_parser.h"
#include "options_info.h"

class options_parser_info : public options_parser
{

public:

    options_parser_info();
    ~options_parser_info();

    bool parse( int argc, char ***argv, options *opts );

};

#endif // OPTIONS_PARSER_INFO_HH_INCLUDED
