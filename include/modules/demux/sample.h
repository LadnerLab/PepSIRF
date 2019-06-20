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
     * A sample can be identified by the first two entries in a tab-separated line. 
     * The first entry is the first identifier, and the second entry is identified similarly. 
     **/
    std::pair<std::string, std::string> string_ids;

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
            string_ids = std::make_pair( id1, id2 );
            name       = sample_name;
            id         = sample_id;
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
     * Get the first id of a sample.
     * @returns const Reference to a sample's first id
     **/
    const std::string& get_first_id() const
        {
            return std::get<0>( string_ids );
        }

    /**
     * Get the second id of a sample.
     * @returns const Reference to a sample's second id
     **/
    const std::string& get_second_id() const 
        {
            return std::get<1>( string_ids );
        }

    /**
     * Get the first id of a sample.
     * @returns Reference to a sample's first id
     **/
     std::string& get_first_id() 
        {
            return std::get<0>( string_ids );
        }

    /**
     * Get the second id of a sample.
     * @returns Reference to a sample's second id
     **/
     std::string& get_second_id()  
        {
            return std::get<1>( string_ids );
        }

     /**
      * Get the ids of this sample. If either (or both) of this sample's ids are 
      * undefined, then we concatenate empty strings, so only the defined sample ids are 
      * included. 
      * @returns std::string The first id of this sample concatenated with the second.
      **/
     std::string get_ids()
         {
             std::string out_str = get_first_id() + get_second_id();
             return out_str;
         }

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

namespace std
{
template <>
struct hash<sample>
{
    std::size_t operator()( const sample& s ) const
        {
            using std::hash;
            using std::string;

            std::string hash_string = s.get_first_id() + s.get_second_id();

            return ( (hash<string>()( hash_string ) ) );
        }

};
}

#endif // SAMPLE_HH_INCLUDED
