#include "metadata_map.h"

std::string metadata_map::build_map( std::string metadata_fname, std::string sequence_name )
    {
        // verify file
        std::cout << "WARNING: Metadata file has been provided and will be implemented with taxonomic ID index ignored." << std::endl;
        //find and verify existence and order for file name/path, pep seq name, spec name : eg. taxtweak_2019-09-12.metadata,Name,Species
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
        std::unordered_map<std::string, std::string> meta_map;
        std::string line;
        std::vector< std::string > metadata_header_row;
        std::getline( metadata_file, line );
        /*std::vector<std::string>::iterator opt_iter;
        opt_iter = ++metadata_options.begin(); // step past file name
        */
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
        std::size_t name_index;
        std::size_t spec_index;
        std::size_t count = 0;
        for( const auto& column_val : metadata_header_row )
            {
                if( column_val.compare( metadata_options[ 1 ] ) == 0 )
                    name_index = count;
                if( column_val.compare( metadata_options[ 2 ] ) == 0 )
                    spec_index = count;
                count++;
            }
    // store name column and species column info in map
        std::vector< std::string > metadata_row;
        while( std::getline( metadata_file, line) )
            {
                boost::split( metadata_row, line, boost::is_any_of( "\t" ) );
                meta_map.insert( std::make_pair( metadata_row.at( name_index ), metadata_row.at( spec_index ) ) );
            }
    // pass map to metadata value class to handle value retrieval.
        metadata_retrieve mr;
        return mr.find_sequence( meta_map, sequence_name );
    }
