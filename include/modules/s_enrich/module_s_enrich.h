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

class module_s_enrich : public module
{
public:

    void run( options *opts );

};

#endif // MODULE_S_ENRICH_HH_INCLUDED 
