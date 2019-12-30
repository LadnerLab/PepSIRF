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
        double geom_mean( Iterator begin,
                          Iterator end
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

    /**
     * Return the sum of the square differences for every item in 
     * range [begin,end).
     * @param begin The first item in the range
     * @param end the last item in the range
     * @param subtrahend The item to subtract from each value in 
     *        range [begin,end)
     * @returns sum( ( x - subtrahend )^2 for x in [begin,end) )
     **/
    template<typename Iterator>
        double squared_diff( Iterator begin,
                             Iterator end,
                             double subtrahend
                           )
        {
            double accumulate = 0;
            for( auto idx = begin; idx != end; ++idx )
                {
                    double diff = *idx - subtrahend;
                    accumulate += ( diff ) * ( diff );
                }
            return accumulate;
        }

    /**
     * Given the population mean, calculate the population 
     * standard deviation of the items in the range [begin, end).
     * @param begin The first item in the range
     * @param end The last item in the range
     * @param the arithmetic mean of the items in the range [begin, end)
     * @returns the population standard deviation.
     **/
    template<typename Iterator>
        double stdev( Iterator begin,
                      Iterator end,
                      double mean
                    )
        {
            double diff = squared_diff( begin, end, mean );
            std::size_t N = std::distance( begin, end );
            return std::sqrt( ( 1 / (double) N ) * diff );
        }

    /**
     * Calculate the population standard deviation of 
     * the range [begin, end).
     * @param begin The first item in the range.
     * @param end The last item in the range.
     * @returns the population standard deviation of the range [begin,end).
     **/
    template<typename Iterator>
        double stdev( Iterator begin,
                      Iterator end
                    )
        {
            double mean = arith_mean( begin, end );
            return stdev( begin, end, mean );;
        }

    /**
     * Calculate the z-score for a value, given the mean and 
     * standard dataset of the data the value is from.
     **/
    template<typename T>
    double zscore( T value, double mean, double stdev )
    {
        return ( value - mean ) / stdev;
    }

    /**
     * Calculate the z-score for each item in the 
     * range [begin,end), outputting to result.
     * @tparam Iterator input iterator to read data from.
     * @tparam OutputIterator the type of iterator to write
     *         results to.
     * @param mean The population arithmetic mean of 
     *        the items in the input range.
     * @param stdev the population standard deviation of the data
     **/
    template<typename Iterator,
             typename OutputIterator>
        void zscore( Iterator begin,
                     Iterator end,
                     OutputIterator result,
                     double mean,
                     double stdev
                   )
        {
            for( auto idx = begin; idx != end; ++idx )
                {
                    *result = zscore( *idx, mean, stdev );
                    ++result;
                }
        }

    /**
     * Calculate the z-score for each item in the range [begin, end), 
     * outputting to result.
     * @param src_begin The first item in the range.
     * @param src_end The last item in the range
     * @param output The output iterator. 
     **/
    template<typename Iterator, typename OutputIterator>
        void zscore( Iterator src_begin,
                     Iterator src_end,
                     OutputIterator output
                   )
        {
            double mean = arith_mean( src_begin, src_end );
            double stand_dev = stdev( src_begin, src_end, mean );

            zscore( src_begin, src_end, output, mean, stand_dev );
        }

};

#endif // STATS_HH_INCLUDED
