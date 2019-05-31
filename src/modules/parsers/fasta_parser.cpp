#include "fasta_parser.h"

fasta_parser::fasta_parser() = default;

std::vector<sequence> fasta_parser::parse( std::string filename )
{
    const char SEQ_START_TOKEN = '>';

    std::vector<sequence> seq_vector;
    sequence seq;

    std::string line;
    std::ifstream input_file( filename, std::ios_base::in );

    while( std::getline( input_file, line ).good() )
        {
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
                            seq = sequence();
                            seq.name = line;
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


