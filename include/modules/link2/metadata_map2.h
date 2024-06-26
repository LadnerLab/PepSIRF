#ifndef METADATA_MAP2_HH_INCLUDED
#define METADATA_MAP2_HH_INCLUDED
#include <string>
#include <vector>
#include <boost/algorithm/string.hpp>
#include <boost/regex.hpp>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>

#include "logger.h"

class metadata_map2
{
 public:

    /**
     * Constructs an unordered map to be used to find species id.
     * @param meta_map metadata retrieve address that will be used in species ID retrieval.
     * @param metadata_fname provides filename, name, and species. Delimited by commas.
    **/
    void build_map( std::unordered_map<std::string, std::unordered_set<std::string>> *meta_map, std::string metadata_fname );

};

#endif // METADATA_MAP2_HH_INCLUDED
