#ifndef PEPTIDE_SCORE_HH_INCLUDED
#define PEPTIDE_SCORE_HH_INCLUDED

template<typename PepType>
class peptide_score
{
 public:
    PepType peptide;
    double zscore;
    double norm_score;
    double raw_score;
    
    peptide_score() = default;
    peptide_score( const PepType pep,
                   const double zscore,
                   const double norm_score,
                   const double raw_score
                 )
        : peptide{ pep },
          zscore{ zscore },
          norm_score{ norm_score },
          raw_score{ raw_score } {}

    peptide_score( const PepType pep,
                   const double zscore,
                   const double norm_score
                 )
        : peptide_score{ pep, zscore, norm_score, 0 }
           {}





};

#endif // PEPTIDE_SCORE_DATA_HH_INCLUDED
