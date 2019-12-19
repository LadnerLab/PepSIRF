#ifndef STATS_HH_INCLUDED
#define STATS_HH_INCLUDED
#include <cmath>
#include <algorithm>
#include <vector>

/**
 * Includes methods to assist in the statistical
 * analyses performed by PepSIRF
 **/
namespace stats
{
    /**
     * Compute the geometric mean of a set of data of type T.
     * @pre T must be a numeric type (std::size_t, double, int, etc.)
     * @pre all of the values in data are greater than zero
     * @note Only non-zero values are considered in the calculation, so 
     *       the 'size' of data is considered the number of non-zero values in 
     *       data.
     * @param data The dataset whose geometric mean to compute.
     **/
    template <typename Iterator
             >
        double geom_mean( const Iterator& begin,
                          const Iterator& end
                        )
        {
            std::size_t n = 0;
            double log_sum = 0;

            for( auto index = begin; index != end; ++index )
                {
                    log_sum += std::log( *index == 0 ? 1 : *index );
                    n += *index != 0;
                }

            log_sum /= (double) n;

            return n == 0 ? 0 : std::exp( log_sum );
        }
};

#endif // STATS_HH_INCLUDED
