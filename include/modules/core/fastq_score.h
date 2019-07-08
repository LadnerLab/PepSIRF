#include "fastq_sequence.h"

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

};

 double get_avg_score( std::string::iterator start, std::string::iterator end )
 {
     int base = phred::phred_base::PHRED33;
     return get_avg_score( start, end, base );
  }

 double get_avg_score( std::string::iterator start, std::string::iterator end,
                       int base
                     )
 {
     for( auto current = start; current != end; ++current )
         {
             sum += *current - base;
         }
     return sum / scores.length();
 }

}
