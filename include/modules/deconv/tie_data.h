#ifndef TIE_DATA_HH_INCLUDED
#define TIE_DATA_HH_INCLUDED
/**
 * Used for data that is relevant to determining and 
 * communicating data about ties.
 **/
namespace tie_data
{
    /**
     * Enum declaring the types of ties that 
     * can be considered. Here we consider 
     * 3 types.
     **/
    enum tie_type
    {
        /**
         * A single-way tie happens when only one species 
         * meets the criteria that considers ties.
         **/
        SINGLE_WAY_TIE = 0,

        /**
         * A two-way tie happens when exactly two species
         * meet the criteria that considers ties.
         **/
        TWO_WAY_TIE,

        /**
         * A k-way tie happens when at least three species
         * meet the criteria that considers ties.
         **/
        K_WAY_TIE
    };

};

#endif // TIE_DATA_HH_INCLUDED
