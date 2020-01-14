#ifndef PREDICATE_HH_INCLUDED
#define PREDICATE_HH_INCLUDED

namespace predicate
{
    bool conjunction( bool p, bool q );
    bool disjunction( bool p, bool q );
    bool implication( bool p, bool q );
    bool biconditional( bool p, bool q );
}; // namespace predicate

#endif // PREDICATE_HH_INCLUDED
