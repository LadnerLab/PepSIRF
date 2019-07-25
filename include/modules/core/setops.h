#ifndef SETOPS_HH_INCLUDED
#define SETOPS_HH_INCLUDED
#include <vector>
#include "maps.h" // used internally

namespace setops
{

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
        


}; // namespace setops


#endif // SETOPS_HH_INCLUDED
