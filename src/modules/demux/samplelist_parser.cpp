#include "samplelist_parser.h"
#include <sstream>

std::vector<sample> samplelist_parser::parse( const std::string filename )
{
    const int FORWARD_ONLY_SIZE        = 2;
    const int FORWARD_AND_REVERSE_SIZE = 3;

    std::ifstream in_stream( filename );

    std::string id1  = "";
    std::string id2  = "";
    std::string name = "";

    std::string line;
    std::size_t line_no = 0;

    std::vector<std::string> split_line;

    int sample_id = 0;

    std::vector<sample> vec;

    while( std::getline( in_stream, line ) )
        {
            boost::split( split_line, line, boost::is_any_of( "\t" ) );
            ++line_no;

            if( split_line.size() == FORWARD_ONLY_SIZE )
                {
                    id1  = split_line[ 0 ];
                    name = split_line[ 1 ];
                }
            else if( split_line.size() == FORWARD_AND_REVERSE_SIZE )
                {
                    id1  = split_line[ 0 ];
                    id2  = split_line[ 1 ];
                    name = split_line[ 2 ];
                }
            else
                {
                    std::stringstream err_str;
                    err_str << "The samplelist file is not formatted correctly!\n" << 
                        "the error occurs here: \n" << line_no << " " <<  line << "\n";
                    throw std::runtime_error( err_str.str() );
                }

            sample samp( id1, id2, name, sample_id );
            vec.push_back( samp );
            ++sample_id;
        }

    return vec;
}
