#include "metadata_map.h"
#include "algorithm"

void build_map( std::string metadata_fname )
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
        std::unordered_map<std::string,
                                   std::unordered_set<scored_entity<std::string,double>>>&
                                   scores_map;
        std::string row;
        std::getline( metadata_file, row );
        std::vector<std::string>::iterator opt_iter;
        opt_iter = ++metadata_options.begin(); // step past file name
        std::for_each( opt_iter, metadata_options.end(), [&]( std::string header )
                {
                    if( row.find( header ) == std::string::npos )
                        {
                            throw std::runtime_error( "Header '" + header + "' could not be found in metadata file. "
                                                    "Verify names being entered to names in metadata file.\n" );
                        }
                }
            );
        std::size_t name_index = std::count( row[ 0 ], row.find( metadata_options[ 1 ] ), "\t" );
        std::size_t spec_index = std::count( row[ 0 ], row.find( metadata_options[ 2 ] ), "\t" );
    // parse metadata to construct map - each row column of name, first, column of spec, second
        while( std::getline( metadata_file, row ) )
            {
                scores_map.insert( std::make_pair( row.substr( name_index, row.find( "\t", name_index ) ),
                                    row.substr( spec_index, row.find( "\t", spec_index ) ) ) );
            }
    // NOTE: What if no string is given for a row id? Currently, it will be an empty string in map.
        metadata_parser mp;
        mp.find_data( retriever );
    }
