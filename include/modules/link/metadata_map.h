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

class metadata_map
{
 public:

    /**
     * Constructs an unordered map to be used to find species id.
     * @param mr metadata retrieve object that will be used in species ID retrieval.
     * @param metadata_fname provides filename, name, and species. Delimited by commas.
    **/
    void build_map( std::string metadata_fname );
   /**
     * Get the species ID using the unordered map built from the metadata file given.
     * @param metadata_map unordered map created from metadata file in metadata_map class.
     * @param sequence_data line of data provided from module_link that contains name and species ID.
    **/
    std::string get_id( std::string sequence_data );

    /**
     * Sets column index for name to be used as a reference for retrieving data from column rows.
     * @param new_val the value to be set as index.
    **/
    void set_name_index( std::size_t new_val );

    /**
     * Gets column index for name to be used as a reference for retrieving data from column rows.
    **/
    std::size_t get_name_index();

    /**
     * Sets column index for species to be used as a reference for retrieving data from column rows.
     * @param new_val the value to be set as index.
    **/
    void set_spec_index( std::size_t new_val );

    /**
     * Gets column index for species to be used as a reference for retrieving data from column rows.
    **/
    std::size_t get_spec_index();

 private:
    /**
     * metadata map used to hold metadata file data.
    **/
    std::unordered_map<std::string,std::string> meta_map;
    /**
     * Column index for name in metadata file.
    **/
    std::size_t name_index;
    /**
     * Column index for species in metadata file.
    **/
    std::size_t spec_index;
};

#endif // METADATA_MAP_HH_INCLUDED
