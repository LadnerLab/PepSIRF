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

    /**
     * Constructs an unordered map to be used to find species id.
     * @param sequence_data set of data that contains name and species.
     * @param metadata_fname provides filename, name, and species. Delimited by commas.
    **/
    std::string build_map( std::string sequence_data, std::string metadata_fname);

 private:
};

#endif // METADATA_MAP_HH_INCLUDED
