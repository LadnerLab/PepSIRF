#include "samplelist_parser.h"
#include <sstream>
#include <boost/algorithm/string.hpp>

// Pass not just filename, but also the string array of headers
std::vector<sample> samplelist_parser::parse( const options_demux *d_opts, const std::vector<flex_idx> flexible_idx_data )
{ 
    std::ifstream samplelist_stream( d_opts->samplelist_fname );
    std::size_t samplename_idx;
    std::vector<std::size_t> index_cols;
    std::string line, header_row;
    std::vector<std::string> split_line;
    std::unordered_set<std::string> sample_id_set;
    std::unordered_set<std::string> index_id_set;
    std::string name;
    int sample_id = 0;
    bool sname_found = false;
    bool index_found = false;
    std::vector<sample> vec;
    std::map<std::string, std::size_t> names;
    std::map<std::string, std::size_t> id_sets;
    bool duplicate_name = false;
    bool duplicate_id = false;

    if( !samplelist_stream.is_open() )
        {
            throw std::runtime_error( "File could not be opened. Verify sample list file exists." );
        }

    std::getline( samplelist_stream, header_row );
    boost::trim_right( header_row );
    boost::split( split_line, header_row, boost::is_any_of( "\t" ) );
    
    for( std::size_t curr_index = 0; curr_index < flexible_idx_data.size(); ++curr_index )
        {
            for( std::size_t curr_col = 0; curr_col < split_line.size(); ++curr_col )
                {
                    if( split_line.at( curr_col ) == flexible_idx_data[curr_index].idx_name )
                        {
                            index_found = true;
                            index_cols.emplace_back( curr_col );
                        }
                }
        }
    for( std::size_t curr_col = 0; curr_col < split_line.size(); ++curr_col )
        {
            if( split_line.at( curr_col ) == d_opts->samplename )
                {
                    sname_found = true;
                    samplename_idx = curr_col;
                }
        }
    if( !sname_found )
        {
            throw std::runtime_error( "The flag \"--sname\" value \'" + d_opts->samplename + "\' could not be found in the sample sheet. "
                                "Verify the sample sheet contains the specified column header names. See demux \"--help\" flag for further help.\n" );
        }
    if( !index_found )
        {
            throw std::runtime_error( "The provided sample sheet does not contain index names found in either the \"--sindex\" option or the \"--fif\" option "
                                      "(depending on which was provided). "
                                      "Verify the sample sheet contains the correct column header names. See demux \"--help\" flag for further information.\n" );
        }
    if( index_cols.size() != flexible_idx_data.size() )
        {
            throw std::runtime_error( "The provided sample sheet does not contain all of the index names provided by either the \"--sindex\" or the "
                                      "\"--fif\" option. Verify the correct indexes are provided and matching for both the \"--sindex\" option or the \"--fif\" option "
                                      "(depending on which is provided).\n" );
        }


    
    // adjust to fit any number of indexes.
    while( std::getline( samplelist_stream, line ) )
        {
            std::string id_set = "";
            boost::trim_right( line );
            boost::split( split_line, line, boost::is_any_of( "\t" ) );
            std::vector<std::string> seqs;
            std::vector<std::string> index_col_ids;
            name = split_line[ samplename_idx ];
            for( const auto& index : index_cols )
                {
                    id_set += split_line[index];
                    sample_id_set.emplace( split_line[index] );
                    index_col_ids.emplace_back( split_line[index] );
                }
            ++names[name];
            ++id_sets[id_set];
            // store the series of sample headers into the sample obj.
            sample samp( index_col_ids, name, sample_id );
            vec.push_back( samp );
            ++sample_id;
        }
    // check for duplicate sample names
    for( auto member : names )
        {                
            if( member.second > 1 )
                {
                    if( !duplicate_name )
                        {
                            std::cout << "WARNING: The following sequence names appear muptiple times" << std::endl;
                            duplicate_name = true;
                        }
                    std::cout << member.first << " Counts: " << member.second << std::endl;
                }
        }
    // check for duplicate id sets
    for( auto member : id_sets )
        {
            if( member.second > 1 )
                {
                    if(!duplicate_id)
                        {
                            std::cout << "WARNING: The following index pairs appear muptiple times" << std::endl;
                            duplicate_id = true;
                        }
                    std::cout << member.first << " Counts: " << member.second << std::endl;
                }
        }

    if( samplelist_stream.bad() )
        {
            throw std::runtime_error( "Encountered error while reading file. Verify sample list file is in .tsv format." );
        }
    std::ifstream index_stream( d_opts->index_fname );
    if( !index_stream.is_open() )
        {
            throw std::runtime_error( "File could not be opened. Verify fasta file exists." );
        }
    while( std::getline( index_stream, line ) )
        {
            boost::trim_right( line );
            if( line.length() > 0 )
                {
                    if( line[0] == '>' )
                        {
                            index_id_set.emplace( line.substr(1) );
                        }
                }
        }
    check_samples( index_id_set, sample_id_set );

    return vec;
}

void samplelist_parser::check_samples( std::unordered_set<std::string>& index_seq_ids, std::unordered_set<std::string>& sample_ids )
    {
        std::vector<std::string> missing_ids;
        for( auto& sample_id : sample_ids )
            {
                if( index_seq_ids.find( sample_id ) == index_seq_ids.end() )
                    {
                        missing_ids.emplace_back( sample_id );
                    }
            }
        if( !missing_ids.empty() )
            {
                std::cout << "WARNING: The following index/barcode names from '--samplelist' are not present in the '--index' file:\n";
                for( auto& id : missing_ids )
                    {
                        std::cout << id << "\n";
                    }
            }

    }
