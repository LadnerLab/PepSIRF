#ifndef SCORED_PEPTIDE_HH_INCLUDED
#define SCORED_PEPTIDE_HH_INCLUDED
#include "peptide.h"

template<typename score_type>
class scored_peptide : public peptide
{
 public:
    scored_peptide() = default;

    scored_peptide( const std::string& name,
                    const std::string& pep,
                    const score_type start_sc
               )
        : peptide( name, pep ), score( start_sc ) {}

    scored_peptide( const std::string& pep,
                    const score_type start_sc
                  )
        : peptide( pep ), score( start_sc ) {}
    


    score_type get_score()
    {
        return score;
    }

    void set_score( score_type new_sc )
    {
        score = new_sc;
    }

 private: 
    score_type score;

};

#endif // SCORED_PEPTIDE_HH_INCLUDED
