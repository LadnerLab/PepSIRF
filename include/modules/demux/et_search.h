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
     * Determine whether multiple best matches for a 
     * sequence have been found. Multiple best matches 
     * means that a single match candidate has two
     * equally-good reference sequences. 
     * @tparam Iterator the iterator type to search 
     *         between.
     * @tparam Get a method that can be used to get items within 
     *         the range [begin, end). For example, if these 
     *         items are pair, a Get method can be implemented that 
     *         retrieves the first item from each pair.
     * @param begin The first item in the range [begin,end)
     * @param end The last item in the range [begin,end)
     * @param ret A method to retrieve items in the range 
     *            [begin,end). 
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

    /**
     * Reference-dependent demultiplexing. For this mode,
     * a reference is used to facilitate error-tolerant demultiplexing.
     * @tparam ref_dep Boolean SFINAE paramater to decide whether reference-independent
     *         demultiplexing should be used. 
     * @param probe The probe to find a match for. This is likely a read from 
     *        the FASTQ file.
     * @param hamming_tol Hamming distance to tolerate when performing the search.
     * @param start The start (0-based) index of the expected read.
     * @param len The length of the expected read.
     * @returns FrequencyVal (FrequencyMap::iterator) pointing to either the 
     *          found sequence in FrequencyMap, or FrequencyMap.end() if the 
     *          reference sequence was not found.
     **/
    template<bool ref_dep = reference_dependent>
    FrequencyVal find( typename std::enable_if<ref_dep, const sequence>::type&
                       probe,
                       std::size_t hamming_tol,
                       std::size_t start,
                       std::size_t len
                     )
    {
        sequence_indexer::query_result
            query_matches;
        sequence *best_match = nullptr;
        sequence seq_temp;

        std::uint32_t num_matches; 

        std::string substr = probe.seq.substr( start, len );

        // Note: hash( sequence& seq ) = hash( seq.seq )
        auto temp = counts.find( sequence( "", substr ) );

        // first check for an exact match in the expected location
        if( temp != counts.end() )
            {
                return temp;
            }

        if( start > 0 )
            {
                substr = probe.seq.substr( start - 1, len );
                temp = counts.find( sequence( "", substr ) );

                if( temp == counts.end()
                    && start > 1
                  )
                    {
                        substr = probe.seq.substr( start - 2, len );
                        temp = counts.find( sequence( "", substr ) );
                    }

                if( temp != counts.end() )
                    {
                        return temp;
                    }
            }

            // shift one to the right, look for exact match.
            // but only if we have enough substring to search
            if( start + 1 + len <= probe.seq.length() )
                {
                    substr = probe.seq.substr( start + 1, len );
                    temp = counts.find( sequence( "", substr ) );
                }

            // look for a match at the expected coordinates within
            // the number of mismatches that are tolerated
            if( hamming_tol > 0 && temp == counts.end() )
                {
                    substr = probe.seq.substr( start, len );
                    seq_temp = sequence( "", substr );
                    num_matches = idx.query( query_matches,
                                             seq_temp,
                                             hamming_tol
                                           );
                    if( num_matches
                        && multiple_best_matches( query_matches.begin(),
                                                  query_matches.end(),
                                                  []( const typename
                                                      sequence_indexer
                                                      ::query_result
                                                      ::value_type& a
                                                    )
                                                  -> std::size_t
                                                  {
                                                      return a.second;
                                                  }
                                                )
                      )
                        {
                            best_match = query_matches[ 0 ].first;
                            temp = counts.find( *best_match );
                        }
                    else if( num_matches )
                        {
                            temp = counts.find( *(query_matches[ 0 ].first) );
                        }
                }

            return temp;
    }

template<typename T>
T *ptr( T& obj )
{
    return &obj;
}

template<typename T>
T *ptr( T *obj )
{
    return obj;
}

/**
 * SFINAE function to enable the return of a 
 * 'indirectionable' object. This means that 
 * performing indirection on the returned value will 
 * return a value of type ptr_maybe. 
 * This may be useful when you need to call the 
 * indirection operator and may need either a pointer to a pointer
 * (when ptr_maybe is a pointer), or a plain pointer (when ptr_maybe is not a pointer). 
 * @tparam ptr_maybe A type that may or may not be pointer. 
 *         
 **/
template<typename ptr_maybe>
    ptr_maybe *make_indirectionable( typename std
                              ::enable_if<std
                              ::is_pointer<ptr_maybe>::value,
                              typename
                              std::remove_pointer<ptr_maybe>::type>::type*&
                              arg
                            )
{
    return &arg;
}

/**
 *
 **/
template<typename ptr_maybe>
    ptr_maybe *make_indirectionable( typename std
                             ::enable_if<!std
                             ::is_pointer<ptr_maybe>::value,
                             typename
                             std::remove_pointer<ptr_maybe>::type>::type*&
                             arg
                           )
{
    return arg;
}




    /**
     * Reference-independent demultiplexing. For this mode, no 
     * reference is used to determine which sequence a read maps to. 
     * Because of this, each read can map to any two sequences:
     * (i) a sequence that's already in the map, and (ii) itself.
     * @tparam ref_dep SFINAE parameter that determines whether 
     *         reference-independent demultiplexing should be used.
     * @pre Each entry in the count map must be the same size. 
     * @param probe The probe whose match in the reference should be found. 
     * @param hamming_tol Dummy parameter, not used. 
     * @param start The start (0-based) index of the expected read.
     * @param len The length of the expected read.
     * @note The hamming distance parameter is ignored, but is maintained so 
     *       reference-independent/dependent demux have the same interface.
     * @returns FrequencyVal (FrequencyMap::iterator) pointing to either the 
     *          found sequence in FrequencyMap, or FrequencyMap.end() if the 
     *          reference sequence was not found.
     **/
    template<bool ref_dep = reference_dependent>
    FrequencyVal find( typename std::enable_if<!ref_dep, const sequence>::type&
                       probe,
                       std::size_t hamming_tol,
                       std::size_t start,
                       std::size_t len
                     )
    {
        std::string sub_str = probe.seq.substr( start, len );

        // stupid hack to the interface doesn't change between
        // reference independent/dependent and we do not get warned
        // about unused parameter
        if( hamming_tol > std::numeric_limits<std::size_t>::max() )
            {
                return counts.end();
            }

        auto iter = counts.find( sequence( "", sub_str ) );
        if( iter == counts.end() )
            {
                using mapped_t_ptr_maybe = typename FrequencyMap::mapped_type;
                using mapped_t = typename std
                    ::remove_pointer<mapped_t_ptr_maybe>::type;

                mapped_t *new_val;

                std::allocator<mapped_t>
                    alloc;
                auto value_size = ptr( counts.begin()->second )->size();

                new_val = alloc.allocate( 1 );

                alloc.construct( new_val,
                                 value_size, 0
                               );

                auto x = counts.emplace( sequence( "", sub_str ),
                                         *make_indirectionable
                                         <mapped_t_ptr_maybe>
                                         (new_val)
                                       );
                iter = x.first;
            }
        return iter;
    }



private:
    const sequence_indexer& idx;
    FrequencyMap& counts;

};

#endif // ET_SEARCH_HH_INCLUDED
