#ifndef ET_SEARCH_HH_INCLUDED
#define ET_SEARCH_HH_INCLUDED
#include "sequence_indexer.h"

/**
 * Error-tolerant (et) sequence search mechanism that 
 * allows for the searching of a data structure with 
 * error tolerance. 
 **/
template<typename FrequencyMap,
         bool reference_independent = false
       >
class et_seq_search
{
public:
    using FrequencyVal = typename FrequencyMap::iterator;

    et_seq_search( const sequence_indexer& ref_idx )
        : idx{ ref_idx }, counts{}
          
    {}

    /**
     *
     **/
    FrequencyVal find( const sequence& probe,
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
