#ifndef OVERLAP_DATA_HH_INCLUDED
#define OVERLAP_DATA_HH_INCLUDED
#include "evaluation_strategy.h"

template<typename T>
class overlap_data
{

 public:
    overlap_data() = default;

    overlap_data( evaluation_strategy::tie_eval_strategy strat,
                  T a_to_b, T b_to_a
                  ) : eval_strat( strat ),
                      a_to_b( a_to_b),
                      b_to_a( b_to_a ) {}

 overlap_data( evaluation_strategy::tie_eval_strategy strat ) : eval_strat( strat )
    {}


 overlap_data( T a_to_b, T b_to_a ) : a_to_b( a_to_b ), b_to_a( b_to_a )
    {}


    bool sufficient( evaluation_strategy::tie_eval_strategy strategy,
                     T threshold
                   )
    {
                return ( threshold == T() )
                    || ( a_to_b >= threshold
                         && b_to_a >= threshold
                       );
    }

    bool sufficient( T threshold )
    {
        return sufficient( eval_strat, threshold );
    }

 private:
    T a_to_b;
    T b_to_a;
    evaluation_strategy::tie_eval_strategy eval_strat;
};



#endif // OVERLAP_DATA_HH_INCLUDED
