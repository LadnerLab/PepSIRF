#ifndef SETOPS_HH_INCLUDED
#define SETOPS_HH_INCLUDED
#include <vector>
#include "maps.h" // used internally

namespace setops
{

    template<class I, class K>
        void intersection( I& dest,
                           const K& first,
                           const K& second
                         )
        {
            for( const auto& elem : first )
                {
                    if( second.find( elem ) != second.end() )
                        {
                            dest.insert( elem );
                        }
                }

        }


    template<class I, class K>
        void intersection( I& dest,
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
                            dest.push_back( it );
                        }
                }


        }
        


}; // namespace setops


#endif // SETOPS_HH_INCLUDED
