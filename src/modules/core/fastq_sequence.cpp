#include "fastq_sequence.h"

double fastq_sequence::get_avg_score()
{
    return get_avg_score( scores.begin(), scores.end() );
}

double fastq_sequence::get_avg_score( std::string::iterator start, std::string::iterator end )
{
    double sum = 0;
    for( auto current = start; current != end; ++current )
        {
            sum += *current - base;
        }
    return sum / scores.length();
}
