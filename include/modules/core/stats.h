#ifndef STATS_HH_INCLUDED
#define STATS_HH_INCLUDED
#include <cmath>
#include <algorithm>
#include <vector>
#include <iostream>
#include <numeric>

/**
 * Includes methods to assist in the statistical
 * analyses performed by PepSIRF
 **/
namespace stats
{
    /**
     * A default epsilon value for determining 'closeness' of two values.
     **/
    constexpr double DEFAULT_EPSILON = 0.001;

    /**
     * Determine the relative difference between 
     * two values of type T.
     **/
    template<typename T>
    double relative_difference( T a, T b )
    {
        T a_abs = std::abs( a );
        T b_abs = std::abs( b );

        T d = std::max( a_abs, b_abs );
        return d == 0.0 ? 0.0 : std::abs( a - b ) / d;
    }

    /**
     * Determine whether the relative_difference between 
     * a and b is less than or equal to the tolerance.
     **/
    template<typename T>
    bool is_close( T a, T b, double tol )
    {
        return relative_difference( a, b ) <= tol;
    }
    
    /**
     * Determine whether the relative difference between 
     * too values is less than or equal to the default 
     * epsilon value.
     **/
    template<typename T>
    bool is_close( T a, T b )
    {
        return is_close( a, b, DEFAULT_EPSILON );
    }
    
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
     * Calculate the arithmetic mean of a range [begin,end).
     * @param begin The first item in the range
     * @param end the last item in the range
     * @returns the arithmetic mean of all items in the range.
     **/
    template<typename Iterator>
        double arith_mean( Iterator begin,
                           Iterator end
                         )
        {
            double sum = std::accumulate( begin, end, 0 );
            return sum / std::distance( begin, end );
        }

    /**
     * Compute the median from a range of items.
     * @param begin The first item in the range 
     * @param end the last item in the range
     * @returns The median of the items in the range [begin,end)
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
