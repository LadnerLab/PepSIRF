#ifndef SCORED_PEPTIDE_HH_INCLUDED
#define SCORED_PEPTIDE_HH_INCLUDED
#include "peptide.h"
#include "scored_entity.h"

/**
 * A scored_peptide is a peptide that has been 
 * associated with a score.
 * The score_type of a peptide is the type that is 
 * used to score a peptide. (int, double, etc.)
 **/
template<typename score_type>
class scored_peptide
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
        : data( peptide( name, pep ), start_sc ) {}

    /**
     * Construct a scored_peptide with a peptide sequence,
     * and score, leaving the name default ("")
     * @param start_sc the starting score for the peptide
     **/
    scored_peptide( const std::string& pep,
                    const score_type start_sc
                  )
        : data( peptide( pep ), start_sc ) {}
    
    /**
     * Return the name of this peptide.
     * @returns constant reference to this peptide's name
     **/
    const std::string& get_name() const
    {
        return data.get_key().get_name();
    }

    /**
     * Return constant reference to this 
     * object's sequence.
     **/
    const std::string& get_sequence() const
    {
        return data.get_key().get_sequence();
    }

    /**
     * Return (by value) the peptide's score.
     * It is expected that score_type is relatively 
     * small, resulting in small overhead for copying,
     * @returns this peptide's score.
     **/
    score_type get_score() const
    {
        return data.get_score();
    }

    /**
     * Set the score of the peptide.
     * @param new_sc The peptide's new score.
     **/
    void set_score( score_type new_sc )
    {
        data.set_score( new_sc );
    }

    /**
     * Determine whether this scored_peptide is less than comp.
     * For two scored_peptides a and b, we say a < b iff
     * a.get_score() < b.get_score()
     * @param comp The peptide with which we are comparing this.
     * @returns boolean true if this < comp, false otherwise
     **/
    bool operator<( const scored_peptide& comp ) const
    {
        return get_score() < comp.get_score();
    }

    /**
     * Determine whether this scored_peptide is greater than comp.
     * For two scored_peptides a and b, we say a > b iff
     * a.get_score() > b.get_score()
     * @param comp The peptide with which we are comparing this.
     * @returns boolean true if this > comp, false otherwise
     **/
    bool operator>( const scored_peptide& comp ) const
    {
        return get_score() > comp.get_score();
    }


 private: 
    /**
     * The score of a peptide, defined 
     * by some scoring metric.
     **/
    scored_entity<peptide,score_type> data;

};

#endif // SCORED_PEPTIDE_HH_INCLUDED
