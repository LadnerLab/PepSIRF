#include "samplelist_parser.h"
#include <sstream>
#include <boost/algorithm/string.hpp>

// Pass not just filename, but also the string for headers.
std::vector<sample> samplelist_parser::parse( const std::string filename )
{
    const int FORWARD_ONLY_SIZE        = 2;
    const int FORWARD_AND_REVERSE_SIZE = 3;

    std::ifstream in_stream( filename );
    // parse the header string if not empty, otherwise use default.
    // Instead of these three columns, headers string will be parsed into a set.
    std::string id1  = "";
    std::string id2  = "";
    std::string name = "";

    std::string line;
    std::size_t line_no = 0;

    std::vector<std::string> split_line;

    int sample_id = 0;

    std::vector<sample> vec;

    if( !in_stream.is_open() )
        {
            throw std::runtime_error( "File could not be opened. Verify sample list file exists." );
        }
    while( std::getline( in_stream, line ) )
        {
            boost::trim_right( line );
            boost::split( split_line, line, boost::is_any_of( "\t" ) );
            ++line_no;
            // Now it will be more nondeterministic what number of headers will be used.
            // We must check to verify the size of split_line is equal to the size of the headers set.
            // Read the first line and assume it to be a line of headers. Store each into an array samplefile_headers.
            // for each header from the input file, check it is found in the set of specified headers. If not throw error.
            // Track every index that is to be included and store for each following line those indexes into the sample object.
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
    if( in_stream.bad() )
        {
            throw std::runtime_error( "Encountered error while reading file. Verify sample list file is in .tsv format." );
        }

    return vec;
}
