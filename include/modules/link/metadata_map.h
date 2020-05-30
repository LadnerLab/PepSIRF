#ifndef METADATA_MAP_HH_INCLUDED
#define METADATA_MAP_HH_INCLUDED

#include <string>
#include <vector>
#include <boost/algorithm/string.hpp>
#include <iostream>
#include <fstream>
#include <algorithm>
#include "metadata_retrieve.h"

class metadata_map
{
    public:
    //function call operator override


    std::string build_map( std::string metadata_fname, std::string sequence_name );

    private:
};

#endif // METADATA_MAP_HH_INCLUDED
