#ifndef SPECIES_DATA_HH_INCLUDED
#define SPECIES_DATA_HH_INCLUDED
#include <string>
#include "species_id.h"
#include "scored_peptide.h"

/**
 * A class to take a 'snapshot' of scoring information for 
 * a species. This class tracks score, count, 
 * and the score of the species' highest-scoring peptide. 
 **/
class species_data
{
 public:
    /**
     * Default constructor.
     **/
    species_data() = default;

    /**
     * Argument constructor.
     * @param species The id for the species whose data is being specified.
     * @param score The score for this species.
     * @param count The count of this species.
     * @param highest_scoring_peptide The peptide whose score is greatest for this species
     * @pre score and count cannot be negative
     **/
    species_data( const species_id<std::string>& species,
                  const double score,
                  const double count,
                  const scored_peptide<double>&
                  highest_scoring_peptide
                ) : spec_id( species ),
                    score( score ),
                    count( count ),
                    highest_scoring_peptide( highest_scoring_peptide )
                    {}

    void set_score( const double new_val );

    void set_count( const double new_val );

    scored_peptide<double>& get_highest_scoring_peptide();

    void set_highest_scoring_peptide( const scored_peptide<double>& new_val );

    double get_score();
    double get_count();

    const scored_peptide<double>& get_highest_scoring_peptide() const;
    const std::string& get_id() const;


 private:
    /**
     * The id of this species, should be unique to a species.
     **/
    species_id<std::string> spec_id;

    /**
     * The score of the species, as defined by some
     * scoring mechanism.
     **/
    double score;

    /**
     * The count of the species, as defined by some
     * counting mechanism. 
     **/
    double count;

    /**
     * A scored peptide representing this species' highest-scoring 
     * peptide. 
     **/
    scored_peptide<double> highest_scoring_peptide;

};


#endif // SPECIES_DATA_HH_INCLUDED
