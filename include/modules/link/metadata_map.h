#ifndef METADATA_MAP_HH_INCLUDED
#define METADATA_MAP_HH_INCLUDED
#include <string>
#include <vector>
#include <boost/algorithm/string.hpp>
#include <boost/regex.hpp>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <unordered_map>

#include "logger.h"

class metadata_map
{
 public:

    /**
     * Constructs an unordered map to be used to find species id.
     * @param meta_map metadata retrieve address that will be used in species ID retrieval.
     * @param metadata_fname provides filename, name, and species. Delimited by commas.
    **/
    void build_map( std::unordered_map<std::string, std::string> *meta_map, std::string metadata_fname );
   /**
     * Get the species ID using the unordered map built from the metadata file given.
     * @param meta_map unordered map address created from metadata file in metadata_map class.
     * @param sequence_data line of data provided from module_link that contains name and species ID.
    **/
   static std::string get_id( std::string sequence_data, std::unordered_map<std::string, std::string> *meta_map );

};

#endif // METADATA_MAP_HH_INCLUDED
