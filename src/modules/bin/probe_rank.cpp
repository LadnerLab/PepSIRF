#include "probe_rank.h"
#include <cmath>


double probe_rank::round_to_factor( const double value ) const
{
    if( !rounding_factor ){ return std::round( value ); }

    const int adjustment = std::pow( 10, rounding_factor );
    return std::floor( value * ( adjustment ) + 0.5 ) / adjustment;
}

void probe_rank::rank_probe( const score_type score,
                             const probe_type probe
                           )
{
    double rounded_score = round_to_factor( score );

    auto insert_loc = ranked_probes.insert( std::make_pair( rounded_score,
                                                            std::vector<probe_type>()
                                                          )
                                          );

    insert_loc.first->second.emplace_back( probe );
}
