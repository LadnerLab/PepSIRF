#include <limits>
#include "fastq_parser.h"

fastq_parser::fastq_parser() = default;

bool fastq_parser::parse( std::ifstream& input_file, std::vector<sequence>& seq_vector, std::size_t max_num_records )
{
    sequence seq;
    const int STRINGS_PER_RECORD = 4;
    size_t max = max_num_records;

    unsigned int count = 0;

    // potentially reduce the amount of copying that needs to be done,
    // if capacity > max nothing is done
    seq_vector.reserve( max );

    // if max is 0 we want to read the entire file at once
    max = max == 0 ? std::numeric_limits<size_t>::max() : max;

    std::string strings[ STRINGS_PER_RECORD ];
    while( count < max && input_file.good() )
        {
            std::getline( input_file, strings[ 0 ] );
            std::getline( input_file, strings[ 1 ] );
            std::getline( input_file, strings[ 2 ] );
            std::getline( input_file, strings[ 3 ] );

            if( input_file.good() )
                {
                    seq = sequence( strings[ 0 ], strings[ 1 ] );
                    seq_vector.push_back( seq );
                    ++count;
                }
        }
    return count == max || max_num_records == 0 ? true : false;
}
