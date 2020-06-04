
#include "metadata_retrieve.h"

std::string metadata_retrieve::get_id( std::unordered_map< std::string, std::string > metadata_map, std::string sequence_data )
{
    std::string sequence_id_codes = sequence_data.substr( sequence_data.find( "OXX=" ) + 4 );
    std::vector<std::string> id_vec;
    boost::split( id_vec, sequence_id_codes, boost::is_any_of( "," ) );
    std::unordered_map<std::string, std::string>::const_iterator name = metadata_map.find( id_vec[ name_index ] );
    if( name != metadata_map.end() )
        {
            return name->second;
        }
    else
        {
            return "";
        }
}

void metadata_retrieve::set_name_index( std::size_t new_val )
{
    name_index = new_val;
}

void metadata_retrieve::set_spec_index( std::size_t new_val )
{
    spec_index = new_val;
}

std::size_t metadata_retrieve::get_name_index()
{
    return name_index;
}

std::size_t metadata_retrieve::get_spec_index()
{
    return spec_index;
}
