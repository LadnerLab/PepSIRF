#ifndef ET_SEARCH_HH_INCLUDED
#define ET_SEARCH_HH_INCLUDED
#include "sequence_indexer.h"
#include <type_traits>

/**
 * Error-tolerant (et) sequence search mechanism that 
 * allows for the searching of a data structure with 
 * error tolerance. 
 * @tparam Frequencymap the type of map that is used to 
 *         track the sequence occurrence frequency in each sample.
 *         Should implement an interface compatible with std::unordered_map
 * @tparam reference_dependent Boolean specifying whether reference-dependent 
 *         demultiplexing should be done.
 *         For reference-dependent mode, at most 6 attempts are made to match 
 *         A sequence to a read region to a reference. They are as follows:
 *         1. Perfect match in the expected region.
 *         2. Perfect match 1 to the left
 *         3. Perfect match 2 to the left
 *         4. Perfect match 1 to the right
 *         5. Perfect match 2 to the right
 *         6. Hamming distance based at expected location.
 *         In reference independent mode, no reference is checked, 
 *         an item is either inserted into the frequency map and its iterator 
 *         is returned, or the iterator of the existing iterator is returned.
 *         
 **/
template<typename FrequencyMap,
         bool reference_dependent = true
       >
class et_seq_search
{
public:
    using FrequencyVal = typename FrequencyMap::iterator;

    /**
     * Initialize the object. 
     * @param ref_idx reference to ref_idx object. 
     *        This indexer should contain at least 
     *        the sequences in the FrequencyMap.
     * @param A map that associates sequences with an ordered 
     *        container of counts for that sequence, one per 
     *        sample.
     * @note To avoid expensive copying of large data structures,
     *       the et_seq_seqrch only maintains references, and does 
     *       NOT maintain the lifetime of the data structures it uses.
     *       It is important to ensure this object does not try to access 
     *       any deleted objects.
     **/
    et_seq_search( const sequence_indexer& ref_idx,
                   FrequencyMap& counts
                 )
        : idx{ ref_idx }, counts{ counts }
          
    {}

    /**
     *
     **/
    template<typename Iterator,
             typename Get
           > 
    bool multiple_best_matches( Iterator begin,
                                Iterator end,
                                Get ret
                              )
    {

        return std::distance( begin, end ) > 1 &&
            ret( *begin ) == ret( *(++begin) );
    }

                       std::size_t hamming_tol,
                       std::size_t start,
                       std::size_t len
                       )
    {

    }


private:
    const sequence_indexer& idx;
    FrequencyMap counts;

};

#endif // ET_SEARCH_HH_INCLUDED
