#ifndef DISTANCE_TOOLS_HH_INCLUDED
#define DISTANCE_TOOLS_HH_INCLUDED

#include "distance_matrix.h"

namespace distance_tools
{
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
    template<template<typename ...> class ContainerType,
        class InputIt,
        class T,
        class DistanceType,
        class DistanceCalc
        >
        void
        all_distances( ContainerType<DistanceType>& dest,
                       InputIt begin, InputIt end,
                       T& compare,
                       DistanceCalc d
                     )
        {
            for( auto it = begin; it < end; ++it )
                {
                    dest.push_back( d( compare, *it ) );
                }
        }

    template<
        class InputIt,
        class DistanceType,
        class DistanceCalc
        >
        void 
        pairwise_distances_symmetry_optimized( distance_matrix<DistanceType>& dest,
                                               InputIt begin,
                                               InputIt end,
                                               DistanceCalc distance
                                             )
            {
                std::size_t index = 0;
                for( auto current = begin; current != end; ++current, ++index )
                    {
                        all_distances( dest[ index ],
                                       begin + index + 1,
                                       end,
                                       *current,
                                       distance
                                     );
                    }
            }

    template<
        class InputIt,
        class DistanceType,
        class DistanceCalc
        >
        void 
        pairwise_distances( distance_matrix<DistanceType>& dest,
                            InputIt begin,
                            InputIt end,
                            DistanceCalc distance
                          )
        {
            std::size_t index = 0;
            for( auto current = begin; current != end; ++current, ++index )
                {
                    all_distances( dest[ index ],
                                   begin,
                                   end,
                                   *current,
                                   distance
                                 );
                }
        }


}; // namespace distance_tools 

#endif // DISTANCE_TOOLS_HH_INCLUDED
