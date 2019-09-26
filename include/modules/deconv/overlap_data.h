#ifndef OVERLAP_DATA_HH_INCLUDED
#define OVERLAP_DATA_HH_INCLUDED

/**
 * Class for storing overlap data. Stores the overlap of 
 * a with respect to b, and the overlap of b with respect to a.
 **/
template<typename T>
class overlap_data
{

 public:
    /**
     * Default constructor.
     **/
    overlap_data() = default;

    /**
     * Argument constructor, takes values for a's overlap with 
     * b and b's overlap with a.
     **/
    overlap_data( T a_to_b, T b_to_a ) : a_to_b( a_to_b),
                                         b_to_a( b_to_a ) {}

    /**
     * Constructor for reflexive overlap.
     * @param a_to_b The overlap of a with respect to b,
     *        which in this instance is also recorded as the 
     *        overlap of b with respect to a. 
     **/
    overlap_data( T a_to_b )
        {
            a_to_b = a_to_b;
            b_to_a = a_to_b;
        }

    /**
     * Determines whether 
     * the overlap is sufficient, 
     * i.e. both a and b have more than 
     * threshold overlap with eachother.
     * @param threshold The threshold at which 
     *        overlap will be considered sufficient.
     **/
    bool sufficient( T threshold )
    {
                return ( threshold == T() )
                    || ( a_to_b >= threshold
                         && b_to_a >= threshold
                       );
    }

    /**
     * Get the overlap of a with respect to b.
     **/
    T get_a_to_b()
    {
        return a_to_b;
    }

    /**
     * Get the overlap of b with respect to a.
     **/
    T get_b_to_a()
    {
        return b_to_a;
    }


 private:
    T a_to_b;
    T b_to_a;
};



#endif // OVERLAP_DATA_HH_INCLUDED
