#ifndef MODULE_S_ENRICH_HH_INCLUDED
#define MODULE_S_ENRICH_HH_INCLUDED
#include "module.h"
#include "options_s_enrich.h"

/**
 * Base case for value_constrained_by
 **/
template<typename Val, typename UnaryPred>
bool value_constrained_by( const Val& arg, const UnaryPred funct )
{
    return funct( arg );
}

/**
 * Template to determine whether a value is constrained by 
 * each in a list of unary predicates. 
 * @tparam Val the type of the value to test
 * @tparam UnaryPred the type of the current unary predicate to test, must 
 *         be a function compatible with  the signature 
 *         []( const Val val ) -> bool
 * @tparam preds the rest of the predicates to test
 * @returns true if every predicate returns true given the value arg
 **/
template<typename Val, typename UnaryPred, typename... Predicates>
bool value_constrained_by( const Val& arg,
                           const UnaryPred funct,
                           Predicates... preds
                         )
{
    constexpr const bool is_same = std::is_same<decltype( funct( arg ) ), bool>::value;

    static_assert( is_same, "funct must be a functor that takes value of type Val "
                            "and returns bool."
                 );

    return funct( arg ) && value_constrained_by( arg, preds... );
}

/**
 * Capture the values in the range [begin, end) for which 
 * pred returns true.
 * @tparam InputIterator The type of the iterator used as input
 * @tparam OutputIterator the type of the iterator used as output
 * @tparam Pred the predicate used to perform the test for each value
 * @pre Pred has a signature 
 *      equivalent to pred( InputIterator::value_type )->bool
 * @param begin The first item in the range [begin,end)
 * @param end The last item in the range [begin,end)
 * @param dest The location to store the values for which pred is true,
 * @note Behavoir undefined if dest is not large enoughto hold all values for which 
 *       pred is true.
 * @returns OutputIterator pointing after the last item that was evaluated 
 *          as true 
 **/
template<typename InputIterator,
         typename OutputIterator,
         typename Pred
         >
OutputIterator valid_for( const InputIterator begin,
                          const InputIterator end,
                          OutputIterator dest,
                          const Pred pred
                        )
{
    using iter_val = typename InputIterator::value_type;
    using pred_ret_type = decltype( pred( std::declval<iter_val>() ) );

    static_assert( std::is_same<pred_ret_type,
                                bool>::value,
                   "Pred must be a unary predicate that returns bool."
                 );

    for( auto curr = begin; curr != end; ++curr )
        {
            if( pred( *curr ) )
                {
                    *dest = *curr;
                    ++dest;
                }
            
        }
    return dest;
}

class module_s_enrich : public module
{
public:

    void run( options *opts );

};

#endif // MODULE_S_ENRICH_HH_INCLUDED 
