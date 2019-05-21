#ifndef OPTIONS_PARSER_FACTORY_HH_INCLUDED 
#define OPTIONS_PARSER_FACTORY_HH_INCLUDED 

#include <boost/program_options.hpp>
#include <algorithm>
#include <string>
#include <stdexcept>

#include "options_parser.h"
#include "options_parser_demux.h"

class options_parser_factory
{
 public:
    options_parser_factory();
    options_parser *create( int argc, char ***argv );
};

#endif // OPTIONS_FACTORY_HH_INCLUDED 
