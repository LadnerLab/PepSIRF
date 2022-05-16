#include "fif_parser.h"
#include <sstream>
#include <boost/algorithm/string.hpp>

std::vector<flex_idx> fif_parser::parse( const std::string fif_fname )
    {
        std::ifstream in_stream( fif_fname );
        std::string line;
        std::vector<std::string> split_line;
        std::vector<flex_idx> ret_idx_data;
        if( !in_stream.is_open() )
            {
                throw std::runtime_error( "File could not be opened. Verify flexible input file (fif) exists.\n" );
            }

        while( std::getline( in_stream, line ) )
            {
                boost::split( split_line, line, boost::is_any_of( "\t" ) );

                if( split_line.size() < 5 )
                    {
                        throw std::runtime_error( "Verify flexible input file contains 5 columns in every row.\n" );
                    }
                if( std::isupper( split_line[1][0] ) )
                    {
                        split_line[1][0] = std::tolower( split_line[1][0] );
                    }
                flex_idx idx_data_elem( split_line[0], split_line[1], split_line[2], split_line[3], split_line[4] );
                ret_idx_data.emplace_back( idx_data_elem );

            }

        return ret_idx_data;
        
    }
