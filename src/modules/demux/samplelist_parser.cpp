#include "samplelist_parser.h"

std::vector<sample> samplelist_parser::parse( const std::string filename )
{
    std::ifstream in_stream( filename );

    std::string id1;
    std::string id2;
    std::string name;

    int sample_id = 1;

    std::vector<sample> vec;

    while( in_stream >> id1 >> id2 >> name )
        {
            sample samp( id1, id2, name, sample_id );
            vec.push_back( samp );
            ++sample_id;
        }

    return vec;
}
