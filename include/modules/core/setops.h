#ifndef SETOPS_HH_INCLUDED
#define SETOPS_HH_INCLUDED
#include <vector>
#include <unordered_set>
#include "maps.h" // used internally

namespace setops
{
    template<typename K,
             typename V>
        struct get_key
        {
            K operator()( const std::pair<K,V>& kv_pair ) const
            {
                return kv_pair.first;
            }
        };

    template<typename K, typename V,
        template<typename...> class Map>
        struct get_value
        {
            K operator()( const typename Map<K,V>::iterator& kv_pair )
            {
                return kv_pair.second;
            }
        };

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
    template<class I, class K, typename Get>
        void set_intersection( I& dest,
                               const K& first,
                               const K& second,
                               const Get retr
                             )
        {
            for( const auto& elem : first )
                {
                    if( second.find( retr( elem ) ) != second.end() )
                        {
                            dest.insert( dest.end(), retr( elem ) );
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
    template<class I, class K, class Get>
        void set_union( I& dest,
                        const K& first,
                        const K& second,
                        Get retr
                      )
    {
        using ValType = typename K::key_type;

        sequential_set<ValType> union_set;

        for( const auto& f : first )
            {
                union_set.insert( retr( f ) );
            }
        for( const auto& s : second )
            {
                union_set.insert( retr( s ) );
            }

        for( const auto& u : union_set )
            {
                dest.insert( u );
            }
    }

    // A u B
    template<class I, class K>
        void set_union( I& dest,
                        const K& first,
                        const K& second
                      )
    {
        std::unordered_set<typename K::value_type> union_set;

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

