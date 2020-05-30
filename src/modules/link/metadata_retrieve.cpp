
#include "metadata_retrieve.h"

std::string metadata_retrieve::find_sequence( std::unordered_map< std::string, std::string > metadata_map, std::string sequence_name )
    {
        std::unordered_map<std::string, std::string>::const_iterator metadata_iter = metadata_map.find( sequence_name );
        if( metadata_iter != metadata_map.end() )
            {
                std::cout << metadata_map.find( sequence_name )->first << metadata_map.find( sequence_name )->second << std::endl;
                return metadata_map.find( sequence_name )->second;
            }
        else
            {
                std::cout << "Sequence name \"" << sequence_name << "\"  could not be found in metadata file." << std::endl;
                return "";
            }

    }
