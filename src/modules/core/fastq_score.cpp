#include "fastq_score.h"

double fastq_score::get_avg_score( std::string::iterator start,
                                   std::string::iterator end,
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
