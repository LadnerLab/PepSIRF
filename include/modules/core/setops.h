#ifndef SETOPS_HH_INCLUDED
#define SETOPS_HH_INCLUDED
#include <vector>
#include "maps.h" // used internally

namespace setops
{

    // A - B
    template<class I, class K>
        void set_intersection( I& dest,
                               const K& first,
                               const K& second
                             )
        {
            for( const auto& elem : first )
                {
                    if( second.find( elem ) != second.end() )
                        {
                            dest.insert( dest.end(), elem );
                        }
                }

        }

    // A - B
    template<class I, class K>
        void set_intersection( I& dest,
                               const std::vector<K>& first,
                               const std::vector<K>& second
                             )
        {
            sequential_set<K> intersection;

            for( const auto& it : first )
                {
                    intersection.emplace( it );
                }

            for( const auto& it : second )
                {
                    if( intersection.find( it ) != intersection.end() )
                        {
                            dest.insert( dest.end(), it );
                        }
                }


        }

    // A u B
    template<class I, class K>
        void set_union( I& dest,
                        const K& first,
                        const K& second
                      )
    {
        sequential_set<K> union_set;

        for( const auto& f : first )
            {
                union_set.insert( f );
            }
        for( const auto& s : second )
            {
                union_set.insert( s );
            }

        for( const auto& u : union_set )
            {
                dest.insert( u );
            }
    }

    // A - B
    template<class I, class K>
        void set_difference( I& dest,
                             const K& first,
                             const K& second
                           )
        {
            for( const auto& it : first )
                {
                    if( second.find( it ) == second.end() )
                        {
                            dest.insert( it );
                        }
                }
        }

    // ( A - B ) u ( B - A )
    template<class I, class K>
        void set_symmetric_difference(
                                      I& dest,
                                      const K& first,
                                      const K& second
                                     )
    {
        sequential_set<K> a_minus_b;
        sequential_set<K> b_minus_a;

        set_difference<I,K>( a_minus_b,
                             first,
                             second
                           );
        set_difference<I,K>( b_minus_a,
                             second,
                             first
                           );

        set_union<I,K>( dest,
                        a_minus_b,
                        b_minus_a
                      );
    }


}; // namespace setops


#endif // SETOPS_HH_INCLUDED
