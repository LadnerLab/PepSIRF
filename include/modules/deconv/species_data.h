#ifndef SPECIES_DATA_HH_INCLUDED
#define SPECIES_DATA_HH_INCLUDED
#include <string>
#include "species_id.h"

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

    species_data( const species_id<std::string>& species,
                  const double score,
                  const double count,
                  const double highest_scoring_peptide
                ) : spec_id( species ),
                    score( score ),
                    count( count ),
                    highest_scoring_peptide( highest_scoring_peptide )
                    {}

 private:
    species_id<std::string> spec_id;
    double score;
    double count;
    double highest_scoring_peptide;
};


#endif // SPECIES_DATA_HH_INCLUDED
