#ifndef STATS_HH_INCLUDED
#define STATS_HH_INCLUDED
#include <cmath>
#include <algorithm>
#include <vector>
#include <iostream>

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


    /**
     * Compute the median from a range of items.
     * @param begin The first item in the range 
     * @param end the last item in the range
     * @returns The median of the items in the range [begin, end]
     **/
    template<typename Iterator>
        double median( Iterator begin,
                       Iterator end
                     )
    {
        std::size_t size = std::distance( begin, end );

        size_t n = size / 2;
        std::nth_element(begin, begin + n, end );

        auto vn = *( begin + n );

        if( size % 2 )
            {
                return vn;
            }
        else
            {
                std::nth_element( begin,
                                  begin + n - 1,
                                  end
                                );
                auto vn_min_1 = *( begin + n - 1 );
                return 0.5 * ( vn + vn_min_1 );
            }
    }


};

#endif // STATS_HH_INCLUDED
