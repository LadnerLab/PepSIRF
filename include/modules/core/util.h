#ifndef UTIL_HH_INCLUDED
#define UTIL_HH_INCLUDED

#include <cmath>

namespace util
{
    template<class M, class InputIt>
        void pairs_to_map( M& dest,
                           InputIt first,
                           InputIt last
                         )
    {
        for( auto iter = first; iter != last; ++iter )
            {
                dest.emplace( iter->first, iter->second );
            }

    }

    template<class T>
        T divide( T a, T b )
        {
            assert( b > 0 );
            return a / b;
        }

    template<class T>
        T add( T a, T b )
        {
            return a + b;
        }

    template<class T>
        T subtract( T a, T b )
        {
            return a - b;
        }

    template<class T>
        T multiply( T a, T b )
        {
            return a * b;
        }

    template< class N>
    bool is_integer( N check )
    {
        return std::floor( check )
               == std::ceil( check );
    }

    template<typename T, typename U>
        std::ostream& operator<<( std::ostream& out, const std::pair<T, U>& p )
    {
        out << "[" << p.first << ", " << p.second << "]";
        return out;
    }

    template<template<typename, typename...> class ContainerType,
             typename ValueType,
             typename ... Args>
        void print_container( std::ostream& stream, const ContainerType<ValueType, Args...>& structure )
        {
            for( const auto& item : structure )
                {
                    stream << item << " ";
                }
            stream << "\n";
        }

    /**
     * Calculate the distance between a member
     * of the type T with all items between
     * begin and end. Store the distances in a 
     * structure that holds DistanceType objects.
     * @param A structure that holds members of type DistanceType,
     *        which is the type returned by the Distance
     *        function d. Results will be stored in 
     *        the same order in which they were found.
     * @param begin The first item to compare with.
     * @param end The last item to compare. 
     * @param compare The item whose distance 
     *        from all items between begin and end
     *        (inclusive) is to be computed.
     * @param d The distance func, which takes two 
     *        members of type T and returns their distance,
     *        a member of type DistanceType.
     **/
    template<template<typename...> class DestType,
        class InputIt, class T,
        class Distance,
        class DistanceType
        >
        void
        all_distances( DestType<DistanceType>& dest,
                       InputIt begin, InputIt end,
                       T& compare,
                       Distance d
                     )
        {
            for( auto it = begin; it < end; ++it )
                {
                    dest.emplace_back( d( compare, *it ) );
                }
        }
};



#endif // UTIL_HH_INCLUDED
