#include "samplelist_parser.h"
#include <sstream>
#include <boost/algorithm/string.hpp>

// Pass not just filename, but also the string array of headers
std::vector<sample> samplelist_parser::parse( const options_demux *d_opts )
{ 
    std::ifstream in_stream( d_opts->samplelist_fname );
    std::size_t line_no = 0;
    std::size_t samplename_idx, id1_col, id2_col;
    std::string line, header_row;
    std::vector<std::string> split_line;
    std::string name, id1;
    std::string id2 = "";
    int sample_id = 0;
    bool sname_found = false;
    bool id1_found = false;
    bool id2_found = false;
    std::vector<sample> vec;

    if( !in_stream.is_open() )
        {
            throw std::runtime_error( "File could not be opened. Verify sample list file exists." );
        }

    std::getline( in_stream, header_row );
    boost::trim_right( header_row );
    boost::split( split_line, header_row, boost::is_any_of( "\t" ) );

    for( std::size_t curr_col = 0; curr_col < split_line.size(); ++curr_col )
        {
        
            if( split_line.at( curr_col ) == d_opts->samplename  )
                {
                    sname_found = true;
                    samplename_idx = curr_col;
                }
            else if( split_line.at( curr_col ) == d_opts->sample_idx1 )
                {
                    id1_found = true;
                    id1_col = curr_col;
                }
            else if( split_line.at( curr_col ) == d_opts->sample_idx2 )
                {
                    id2_found = true;
                    id2_col = curr_col;
                }
                
        }
    if( !sname_found )
        {
            throw std::runtime_error( "Error: The flag \"--sname\" value \'" + d_opts->samplename + "\' could not be found in the sample sheet. "
                                "Verify the sample sheet contains the specified column header names. See demux \"--help\" flag for further help.\n" );
        }
    if( !id1_found )
        {
            throw std::runtime_error( "Error: The flag \"--sindex1\" value \'" + d_opts->sample_idx1 + "\' could not be found in the sample sheet. "
                                "Verify the sample sheet contains the specified column header names. See demux \"--help\" flag for further help.\n" );
        }
    if( !id2_found )
        {
            std::cout << "Warning: The flag \"--sindex2\" value \'" << d_opts->sample_idx2 << "\' could not be found in the sample sheet therefore "
                         "will not be used. By default this optional second index flag is set to \'Index2\'. See the demux \"--help\" flag for "
                         "further information.\n";
        }
    
    while( std::getline( in_stream, line ) )
        {
            boost::trim_right( line );
            boost::split( split_line, line, boost::is_any_of( "\t" ) );
            ++line_no;

            if( !id2_found )
                {
                    name = split_line[ samplename_idx ];
                    id1 = split_line[ id1_col ];
                }
            else
                {
                    name = split_line[ samplename_idx ];
                    id1 = split_line[ id1_col ];
                    id2 = split_line[ id2_col ];
                }
            
            // store the series of sample headers into the sample obj.
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
