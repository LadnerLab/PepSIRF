#ifndef PREDICATE_HH_INCLUDED
#define PREDICATE_HH_INCLUDED
#include <type_traits>

namespace predicate
{
    bool conjunction( bool p, bool q );
    bool disjunction( bool p, bool q );
    bool implication( bool p, bool q );
    bool biconditional( bool p, bool q );

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
     * pred returns true, passing the return value of get to pred.
     * @tparam InputIterator The type of the iterator used as input
     * @tparam OutputIterator the type of the iterator used as output
     * @tparam Pred the predicate used to perform the test for each value
     * @tparam Get A function that taks a value of type InputIterator::value_type
     *         and returns a member of a type for which pred is defined.
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
        typename Pred,
        typename Get
        >
        OutputIterator valid_for( const InputIterator begin,
                                  const InputIterator end,
                                  OutputIterator dest,
                                  const Pred pred,
                                  Get get
                                  )
        {
            using iter_val = typename InputIterator::value_type;
            using get_ret_type = decltype( get( std::declval<iter_val>() ) );
            using pred_ret_type = decltype( pred( std::declval<get_ret_type>() ) );

            static_assert( std::is_same<pred_ret_type,
                           bool>::value,
                           "Pred must be a unary predicate that returns bool."
                           );

            for( auto curr = begin; curr != end; ++curr )
                {
                    if( pred( get( *curr ) ) )
                        {
                            *dest = *curr;
                            ++dest;
                        }
            
                }
            return dest;
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
            auto get = []( const iter_val& val ) -> const iter_val& { return val; };
            return valid_for( begin, end, dest, pred, get );
        }



}; // namespace predicate

#endif // PREDICATE_HH_INCLUDED
