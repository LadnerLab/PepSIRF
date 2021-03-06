#include "species_data.h"
void species_data::set_score( const double new_val )
{
    score = new_val;
}

void species_data::set_count( const double new_val )
{
    count = new_val;
}

scored_peptide<double>& species_data::get_highest_scoring_peptide()
{
    return highest_scoring_peptide;
}

const std::string& species_data::get_id() const
{
    return spec_id.get_id();
}

void species_data
     ::set_highest_scoring_peptide( const scored_peptide<double>& new_val )
{
    highest_scoring_peptide = new_val;
}

double species_data::get_score()
{
    return score;
}

double species_data::get_count()
{
    return count;
}

const scored_peptide<double>&
species_data::get_highest_scoring_peptide() const
{
    return highest_scoring_peptide;
}
