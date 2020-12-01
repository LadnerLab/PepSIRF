#ifndef PAIRED_SCORE_HH_INCLUDED
#define PAIRED_SCORE_HH_INCLUDED

#include <utility>
#include <string>

class paired_score
{

 public:
    using pair = std::pair<double,double>;

    std::string pep_name;
    // score represents either a norm score or zscore pair
    pair score;
    pair raw_score;

 paired_score( const std::string pep_name,
               const pair& score,
               const pair& raw_score
             )
     : pep_name{ pep_name },
       score{ score },
       raw_score{ raw_score } {}

// this paired score represents data where the peptide may be stored separately (eg. map w/ pep name as key)
paired_score( const pair& score,
              const pair& raw_score
            )
     : score{score},raw_score{raw_score} {}

 paired_score( const std::string& pep_name,
               const pair& score
             )
     : paired_score{ pep_name,
                     score,
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
