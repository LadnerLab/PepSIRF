#include "peptide_scoring.h"
#include <fstream>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/join.hpp>
#include <boost/range/adaptor/transformed.hpp>
#include <boost/filesystem.hpp>
#include <sstream>


void peptide_scoring::parse_peptide_scores( peptide_score_data_sample_major& dest,
                                            std::string ifname
                                          )
{
    std::string line;
    std::ifstream in_file( ifname, std::ios_base::in );

    if( in_file.fail() )
        {
            std::stringstream error_msg;
            error_msg << "Unable to open file '" << ifname << "'";
            throw std::runtime_error( error_msg.str() );
        }

    std::vector<std::string> lines_from_file;
    std::vector<std::string> split_line;
    dest.file_name = boost::filesystem::path( ifname )
                     .filename()
                     .string();

    while( std::getline( in_file, line ) )
        {
            boost::trim_right( line );
            lines_from_file.emplace_back( line );
        }

    std::size_t line_count = lines_from_file.size();
    std::size_t num_peptides = line_count - 1; // -1 because of the header
    std::size_t sample_count = std::count( lines_from_file[ 0 ].begin(),
                                           lines_from_file[ 0 ].end(),
                                           '\t'
                                         );
    dest.pep_names.reserve( num_peptides );
    dest.sample_names.reserve( sample_count );

    std::size_t index = 0;

    for( index = 0; index < num_peptides; ++index )
        {
            dest.pep_names.emplace_back( std::string( "" ) );
        }

    // save the sample names to the samplenames vector
    boost::split( split_line, lines_from_file[ 0 ], boost::is_any_of( "\t" ) );
    std::for_each( split_line.begin() + 1, split_line.end(),
                   [&]( const std::string &item )
                   {
                       dest.sample_names.push_back( item );
                   }
                 );

    dest.scores = labeled_matrix<double, std::string>( num_peptides, sample_count );
    for( index = 1; index < lines_from_file.size(); ++index )
        {
            std::size_t assign_index = index - 1;
            // get the name of the peptide
            const std::string& my_str = lines_from_file[ index ];
            std::size_t pos = my_str.find( '\t' );

            dest.pep_names[ assign_index ] = my_str.substr( 0, pos );

            std::size_t begin_pos = pos + 1;
            
            // for j = 0 -> num_samples grab the number
            for( std::size_t count_index = 0; count_index < sample_count; ++count_index )
                {
                    std::size_t end_pos   = my_str.find( '\t', begin_pos );
                    double val = std::strtod( my_str.substr( begin_pos, end_pos - begin_pos ).c_str(),
                                     nullptr
                                   );

                    dest.scores( assign_index, count_index ) = val;
                    begin_pos = end_pos + 1;
                        
                }
        }

    dest.scores.set_row_labels( dest.pep_names );
    dest.scores.set_col_labels( dest.sample_names );
    dest.scores = dest.scores.transpose();

}
    void peptide_scoring::write_peptide_scores( std::string dest_fname,
                                                peptide_score_data_sample_major& data
                                                )
{
    std::ofstream out_file( dest_fname, std::ios_base::out );
    write_peptide_scores( out_file, data ); 
}

void peptide_scoring::write_peptide_scores( std::ostream& output,
                                            peptide_score_data_sample_major& data
                                          )
{
    data.scores = data.scores.transpose();

    output << "Sequence name\t";

    output << boost::algorithm::join( data.sample_names, "\t" ) << "\n";

    for( std::size_t index = 0; index < data.scores.nrows(); ++index )
        {
            output << data.pep_names[ index ] << "\t";
                        
            std::size_t inner_index = 0;

            for( inner_index = 0; inner_index < data.sample_names.size(); ++inner_index )
                {
                    output << data.scores( index, inner_index );

                    if( inner_index < data.sample_names.size() - 1 )
                        {
                            output << "\t";
                        }
                }
            output << "\n";
        }
}


