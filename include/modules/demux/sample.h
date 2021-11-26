#ifndef SAMPLE_HH_INCLUDED
#define SAMPLE_HH_INCLUDED
#include <string>
#include <utility>

/**
 * Data class to represent a sample. 
 **/
class sample
{
 public:
    /**
     * A sample can be identified by the barcode/tag ids in the tab-separated line. 
     * The first entry is the first identifier, and the second entry is identified similarly, and so on. 
     **/
    std::vector<std::string> string_ids;

    /**
     * The sequences for a sample.
     * Positionally, these items correspond with those in 
     * sample::strings_ids.
     * @note The vector can contain only one entry, indicating that 
     *       either only one index is used for this run or 
     *       that sequences are already indexed on R indices.
     **/
    std::vector<std::string> sequences;
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
    sample( std::vector<std::string> index_ids, std::string sample_name, int sample_id )
        {
            string_ids = index_ids;
            name       = sample_name;
            id         = sample_id;
        }

    /**
     * Construct a sample given both ids and sequences, a sample name and id.
     **/
    sample( std::vector<std::string> index_ids, std::vector<std::string> index_sequences, std::string sample_name, int sample_id )
        {
            string_ids = index_ids;
            name       = sample_name;
            id         = sample_id;
            sequences  = index_sequences;
        }

    /**
     * Copy constructor.
     * @param s Sample object to copy.
     **/
    sample( const sample& s )
        {
            string_ids = s.string_ids;
            name       = s.name;
            id = s.id;
        }

    /**
     * Default constructor
     **/
    sample() = default;

     /**
      * Equality operator for samples. 
      * For two samples a and b, we say a == b iff
      * ( a.get_first_id() == b.get_first_id() ) 
      *  and ( a.get_second_id() == b.get_second_id() ) 
      *
      * @param other sample against which to compare sequences
      **/
    bool operator==( const sample& other ) const
    {
        return string_ids == other.string_ids;
    }
};


template <>
struct std::hash<sample>
{
    std::size_t operator()( const sample& s ) const
        {
            using std::hash;
            using std::string;

            std::string hash_string = "";
            for( const auto& id : s.string_ids )
                {
                    hash_string.append( id );
                }

            return ( (hash<string>()( hash_string ) ) );
        }

};

#endif // SAMPLE_HH_INCLUDED
