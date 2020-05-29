#ifndef METADATA_MAP_HH_INCLUDED
#define METADATA_MAP_HH_INCLUDED

#include <string>
#include <boost/algorithm/string.hpp>
#include <iostream>
#include <fstream>
#include <algorithm>
#include "module_link.h"
#include "metadata_parser.h"

class metadata_map
{
    public:
    //function call operator override

    void metadata_map::build_map( std::string metadata_fname, std::string value );

    private:
};

#endif // METADATA_MAP_HH_INCLUDED
