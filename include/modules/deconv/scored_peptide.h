#ifndef SCORED_PEPTIDE_HH_INCLUDED
#define SCORED_PEPTIDE_HH_INCLUDED
#include "peptide.h"

/**
 * A scored_peptide is a peptide that has been 
 * associated with a score.
 * The score_type of a peptide is the type that is 
 * used to score a peptide. (int, double, etc.)
 **/
template<typename score_type>
class scored_peptide : public peptide
{
 public:

    /**
     * Default constructor.
     **/
    scored_peptide() = default;

    /**
     * Construct a scored_peptide with a name, peptide 
     * sequence, and starting score. name and pep are 
     * defined by the 'peptide' class.
     * @param start_sc The starting score of the peptide.
     *        This value will be copied.
     **/
    scored_peptide( const std::string& name,
                    const std::string& pep,
                    const score_type start_sc
               )
        : peptide( name, pep ), score( start_sc ) {}

    /**
     * Construct a scored_peptide with a peptide sequence,
     * and score, leaving the name default ("")
     * @param start_sc the starting score for the peptide
     **/
    scored_peptide( const std::string& pep,
                    const score_type start_sc
                  )
        : peptide( pep ), score( start_sc ) {}
    


    /**
     * Return (by value) the peptide's score.
     * It is expected that score_type is relatively 
     * small, resulting in small overhead for copying,
     * @returns this peptide's score.
     **/
    score_type get_score()
    {
        return score;
    }

    /**
     * Set the score of the peptide.
     * @param new_sc The peptide's new score.
     **/
    void set_score( score_type new_sc )
    {
        score = new_sc;
    }

 private: 
    /**
     * The score of a peptide, defined 
     * by some scoring metric.
     **/
    score_type score;

};

#endif // SCORED_PEPTIDE_HH_INCLUDED
