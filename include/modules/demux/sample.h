#ifndef SAMPLE_HH_INCLUDED
#define SAMPLE_HH_INCLUDED
#include <string>

/**
 * Data class to represent a sample. 
 **/
class sample
{
 public:
    /**
     * The first id of the sample, will be determined by the first entry in a tab-separated line.
     **/
    std::string first_id;

    /**
     * The second id of the sample, will be determined by the second entry in a tab-separated line.
     **/
    std::string second_id;

    /**
     * The third id of the sample, determined by the third entry in a tab-separated line.
     **/
    std::string name;

    /**
     * An integer id of the sample. Note that in a collection of samples the id of each should be unique.
     **/
    int id;

    /**
     * Construct a sample given ids 1 and 2, a sample name and an id for a sample.
     **/
    sample( std::string id1, std::string id2, std::string sample_name, int sample_id )
        {
            first_id  = id1;
            second_id = id2;
            name      = name;
            id        = sample_id;
        }
};

#endif // SAMPLE_HH_INCLUDED
