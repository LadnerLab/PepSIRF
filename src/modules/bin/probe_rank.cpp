#include "probe_rank.h"
#include <cmath>
#include <iostream>

double probe_rank::round_to_factor( const double value ) const
{
    const uint adjustment = std::pow( 10, rounding_factor );
    const double adjusted_val = value * adjustment;

    // round half to even
    if( adjusted_val + 0.5 == std::ceil( adjusted_val ) )
        {
            if( !std::fmod( adjusted_val, 2.0 ) )
                {
                    return adjusted_val / adjustment;
                }
            return ( adjusted_val + 0.5 ) / adjustment;
        }

    return std::ceil( adjusted_val - 0.5 ) / adjustment;
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


probe_rank::rank_track_type::iterator
probe_rank::get_probes_with_rank( const score_type rank )
{
    return ranked_probes.find( round_to_factor( rank ) );
}

probe_rank::rank_track_type::const_iterator
probe_rank::get_probes_with_rank( const score_type rank ) const
{
    return ranked_probes.find( round_to_factor( rank ) );
}


const probe_rank::rank_track_type&
    probe_rank::get_probe_ranks() const
    {
        return ranked_probes;
    }
