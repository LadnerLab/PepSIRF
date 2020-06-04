#include "metadata_map.h"

std::string metadata_map::build_map( std::string sequence_data, std::string metadata_fname)
    {
        std::vector<std::string> metadata_options;
        metadata_options.reserve( 3 );
        boost::split( metadata_options, metadata_fname, boost::is_any_of( "," ) );
        std::ifstream metadata_file( metadata_options[0], std::ios_base::in );
        if( !metadata_file.is_open() )
            {
                throw std::runtime_error( "File could not be opened. Verify metadata file exists.\n" );
            }
        if( metadata_options.size() < 3 )
            {
                throw std::runtime_error( "Missing required specifications for meta flag. Must include metadata file name, "
                                        "sequence name column, and species identification column.\n");
            }
        metadata_retrieve mr;
        std::unordered_map<std::string, std::string> meta_map;
        std::string line;
        std::vector<std::string> metadata_header_row;
        std::getline( metadata_file, line );
        std::for_each( ++metadata_options.begin(), metadata_options.end(), [&]( std::string header )
                {
                    if( line.find( header ) == std::string::npos )
                        {
                            throw std::runtime_error( "Header '" + header + "' could not be found in metadata file. "
                                                    "Verify names being entered to names in metadata file.\n" );
                        }
                }
        );
        boost::split( metadata_header_row, line, boost::is_any_of( "\t" ) );
        std::size_t count = 0;
        for( const auto& column_val : metadata_header_row )
            {
                if( column_val.compare( metadata_options[ 1 ] ) == 0 )
                    mr.set_name_index( count );
                if( column_val.compare( metadata_options[ 2 ] ) == 0 )
                    mr.set_spec_index( count );
                count++;
            }
        std::vector<std::string> metadata_row;
        while( std::getline( metadata_file, line) )
            {
                boost::split( metadata_row, line, boost::is_any_of( "\t" ) );
                meta_map.insert( std::make_pair( metadata_row.at( mr.get_name_index() ), metadata_row.at( mr.get_spec_index() ) ) );
            }
        return mr.get_id( meta_map, sequence_data );
    }
