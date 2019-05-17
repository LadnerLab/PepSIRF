#ifndef OPTIONS_PARSER_HH_INCLUDED
#define OPTIONS_PARSER_HH_INCLUDED
#include <boost/program_options.hpp>
#include <iostream>
#include <iterator>
#include <exception>

#include "options.h"


class options_parser
{
public:
    options_parser();
    void parse( int argc, char ***argv, options& opts );
};





#endif //OPTIONS_PARSER_HH_INCLUDED
