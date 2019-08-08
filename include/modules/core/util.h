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
};



#endif // UTIL_HH_INCLUDED
