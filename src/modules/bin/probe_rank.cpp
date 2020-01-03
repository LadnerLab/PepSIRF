#include "probe_rank.h"
#include <cmath>


double probe_rank::round_to_factor( const double value ) const
{
    if( !rounding_factor ){ return std::round( value ); }

    const int adjustment = std::pow( 10, rounding_factor );
    return std::floor( value * ( adjustment ) + 0.5 ) / adjustment;
}
