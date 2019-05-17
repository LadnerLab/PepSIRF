#ifndef OPTIONS_PARSER_HH_INCLUDED
#define OPTIONS_PARSER_HH_INCLUDED
#include <boost/program_options.hpp>
#include <iostream>
#include <iterator>
#include <exception>


class OptionsParser
{
public:
    OptionsParser();
    void parse( int argc, char ***argv );
};





#endif //OPTIONS_PARSER_HH_INCLUDED
