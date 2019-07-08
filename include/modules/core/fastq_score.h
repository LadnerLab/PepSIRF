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
                     )
 {
     double sum = 0;
     double cnt = 0;
         
     std::for_each( start, end,
                    [&]( const int& curr )
                    {
                        sum += curr - base;
                        ++cnt;
                    }
                  );

     // don't try to divide by 0
     return cnt > 0 ? sum / cnt : 0;
 }

}; //namespace fastq_score

#endif // FASTQ_SCORE_HH_INCLUDED
