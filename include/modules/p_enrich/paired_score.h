#ifndef PAIRED_SCORE_HH_INCLUDED
#define PAIRED_SCORE_HH_INCLUDED

#include <utility>
#include <string>

class paired_score
{

 public:
    using pair = std::pair<double,double>;

    std::string pep_name;
    pair zscore;
    pair norm_score;
    pair raw_score;

 paired_score( const std::string pep_name,
               const pair& zscore,
               const pair& norm_score,
               const pair& raw_score
             )
     : pep_name{ pep_name },
       zscore{ zscore },
       norm_score{ norm_score },
       raw_score{ raw_score } {}

 paired_score( const std::string& pep_name,
               const pair& zscore,
               const pair& norm_score
             )
     : paired_score{ pep_name,
                     zscore,
                     norm_score,
                     { 0.0, 0.0 }
                   } {}

       friend std::ostream& operator<<( std::ostream& stream,
                                        const paired_score& pair
                                      )
       {
           stream << pair.pep_name;
           return stream;
       }

};

#endif // PAIRED_SCORE_HH_INCLUDED
