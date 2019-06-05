#include "samplelist_parser.h"
#include <iostream>

std::vector<sample> samplelist_parser::parse( const std::string filename )
{
    const int FORWARD_ONLY_SIZE        = 2;
    const int FORWARD_AND_REVERSE_SIZE = 3;

    std::ifstream in_stream( filename );

    std::string id1  = "";
    std::string id2  = "";
    std::string name = "";

    std::string line;

    std::vector<std::string> split_line;

    int sample_id = 0;

    std::vector<sample> vec;

    while( std::getline( in_stream, line ) )
        {
            boost::split( split_line, line, boost::is_any_of( "\t" ) );

            if( split_line.size() == FORWARD_ONLY_SIZE )
                {
                    id1  = split_line[ 0 ];
                    name = split_line[ 2 ];
                }
            else if( split_line.size() == FORWARD_AND_REVERSE_SIZE )
                {
                    id1  = split_line[ 0 ];
                    id2  = split_line[ 1 ];
                    name = split_line[ 2 ];
                }
            else
                {
                    throw std::runtime_error( "The samplelist file is not formatted correctly!" );
                }

            sample samp( id1, id2, name, sample_id );
            vec.push_back( samp );
            ++sample_id;
        }

    return vec;
}
