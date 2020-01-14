#include "predicate.h"

namespace predicate
{
    bool biconditional( bool p, bool q )
    {
        return implication( p, q ) && implication( q, p );
    }

    bool implication( bool p, bool q )
    {
        return !p || q;
    }

    bool conjunction( bool p, bool q )
    {
        return p && q;
    }

    bool disjunction( bool p, bool q )
    {
        return p || q;
    }

}; // namespace predicate
