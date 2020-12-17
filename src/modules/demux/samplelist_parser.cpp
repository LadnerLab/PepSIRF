#include "samplelist_parser.h"
#include <sstream>
#include <boost/algorithm/string.hpp>

// Pass not just filename, but also the string array of headers
std::vector<sample> samplelist_parser::parse( const std::string filename, std::unordered_set<std::string> header_names )
                                                                           // A question I need to consider for this process is:
                                                                           // Do I need to restructure how samplelist data is stored?
{                                                                          //  - I know the order will not change per column being used
                                                                           //  - This is a pretty straightforward update. Each usable column is specified and verified to exist.
    std::ifstream in_stream( filename );
    std::size_t line_no = 0;
    std::string line;
    std::vector<std::string> split_line;
    std::vector<std::string> samplelist_data;
    int sample_id = 0;

    std::vector<sample> vec;

    if( !in_stream.is_open() )
        {
            throw std::runtime_error( "File could not be opened. Verify sample list file exists." );
        }

    // Default process to obtaining samples from samplelist
    if( header_names.empty() )
        {
            const int ID1_ONLY_SIZE     = 2;
            const int ID2_INCLUDED_SIZE = 3;
            // Assume 2-3 columns: ( idx1 or ( idx1 and idx2 )) and samplename
            
            while( std::getline( in_stream, line ) )
                {
                    boost::trim_right( line );
                    boost::split( split_line, line, boost::is_any_of( "\t" ) );
                    ++line_no;

                    if( split_line.size() == ID1_ONLY_SIZE )
                        {
                            samplelist_data.emplace_back( split_line[ 0 ] ); // id1
                            samplelist_data.emplace_back( split_line[ 1 ] ); // id2
                        }
                    else if( split_line.size() == ID2_INCLUDED_SIZE )
                        {
                            samplelist_data.emplace_back( split_line[ 0 ] ); // id1
                            samplelist_data.emplace_back( split_line[ 1 ] ); // id2
                            samplelist_data.emplace_back( split_line[ 2 ] ); // samplename
                        }
                    else
                        {
                            std::stringstream err_str;
                            err_str << "The samplelist file is not formatted correctly!\n" <<
                                "the error occurs here: \n" << line_no << " " <<  line << "\n";
                            throw std::runtime_error( err_str.str() );
                        }
                    // store the series of sample headers into the sample obj.
                    sample samp( samplelist_data, sample_id );
                    vec.push_back( samp );
                    ++sample_id;
                }
        }
    else if( !header_names.empty() ) // user-defined inclusivity of headers -> store samples only for columns provided by header_names
        {
            std::string header_row, name;
            std::size_t col_idx;
            std::vector<int> column_idxs;

            std::getline( in_stream, header_row );
            boost::trim_right( header_row );
            boost::split( split_line, header_row, boost::is_any_of( "\t" ) );

            // get list of indexes of names col header names matching --header_names to samplelist col headers 
            for( col_idx = 0; col_idx < split_line.size(); ++col_idx )
                {
                    name = split_line.at( col_idx );
                    if( header_names.find( name ) != header_names.end() )
                        {
                            column_idxs.emplace_back( col_idx );
                        }
                }
            // Reading samplelist file starting at second row
            while( std::getline( in_stream, line ) )
                {
                    boost::trim_right( line );
                    boost::split( split_line, line, boost::is_any_of( "\t" ) );
                    ++line_no;
                    for( const auto& col : column_idxs )
                        {
                            samplelist_data.emplace_back( split_line[ col ] );
                        }
                    sample samp( samplelist_data, sample_id );
                    vec.push_back( samp );
                    ++sample_id;
                }
        }
    if( in_stream.bad() )
        {
            throw std::runtime_error( "Encountered error while reading file. Verify sample list file is in .tsv format." );
        }

    return vec;
}
