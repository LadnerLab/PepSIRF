#include "fasta_parser.h"
#include <boost/algorithm/string.hpp>

fasta_parser::fasta_parser() = default;

std::vector<sequence> fasta_parser::parse( std::string filename )
{
    const char SEQ_START_TOKEN = '>';

    std::vector<sequence> seq_vector;
    sequence seq;

    std::string line;
    std::ifstream input_file( filename, std::ios_base::in );

    if( !input_file.is_open() )
        {
            throw std::runtime_error( "File could not be opened. Verify fasta file exists." );
        }
    while( std::getline( input_file, line ) )
        {
            boost::trim_right( line );
            // skip any empty lines
            if( line.length() > 0 )
                {

                    if( line[ 0 ] == SEQ_START_TOKEN )
                        {
                            if( seq.name.compare( "" ) )
                                {
                                    // if the name has been changed, add the last-created sequence
                                    seq_vector.push_back( seq );
                                }

                            // create a new sequence, set its name to the line
                            #pragma GCC diagnostic push
                            #pragma GCC diagnostic ignored "-Wdeprecated-copy"
                            seq = sequence();
                            seq.name = line.substr( 1 );
                            #pragma GCC diagnostic pop
                        }
                    else
                        {
                            seq.seq.append( line );
                        }
                }

        }

    // last sequence would otherwise not have been added
    seq_vector.push_back( seq );

    return seq_vector;
}
