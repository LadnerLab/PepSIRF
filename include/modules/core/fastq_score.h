#ifndef FASTQ_SCORE_HH_INCLUDED
#define FASTQ_SCORE_HH_INCLUDED

#include "fastq_sequence.h"
#include <algorithm>

namespace fastq_score
{
    
namespace phred
{

    constexpr int PHRED33_ASCII_BASE = 33;
    constexpr int PHRED64_ASCII_BASE = 64;

    enum phred_base
    {
        PHRED33 = PHRED33_ASCII_BASE,
        PHRED64 = PHRED64_ASCII_BASE
    };

}; // namesapce phred

 double get_avg_score( std::string::iterator start, std::string::iterator end,
                       int base
                     );

}; //namespace fastq_score

#endif // FASTQ_SCORE_HH_INCLUDED
