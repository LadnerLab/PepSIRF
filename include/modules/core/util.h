#ifndef UTIL_HH_INCLUDED
#define UTIL_HH_INCLUDED

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
};



#endif // UTIL_HH_INCLUDED
