#ifndef PAIRED_SCORE_HH_INCLUDED
#define PAIRED_SCORE_HH_INCLUDED
#include <utility>

class paired_score
{

 public:
    using pair = std::pair<double,double>;
    pair zscore;
    pair norm_score;
    pair raw_score;

 paired_score( const pair& zscore,
               const pair& norm_score,
               const pair& raw_score
             )
     : zscore{ zscore },
       norm_score{ norm_score },
       raw_score{ raw_score } {}

 paired_score( const pair& zscore,
               const pair& norm_score
             )
     : paired_score{ zscore,
                     norm_score,
                     { 0.0, 0.0 }
                   } {}

};

#endif // PAIRED_SCORE_HH_INCLUDED
