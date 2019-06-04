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
    std::size_t operator()( sample& s ) const
        {
            using std::hash;
            using std::string;

            std::string hash_string = s.get_first_id() + s.get_second_id();

            return ( (hash<string>()( hash_string ) ) );
        }

};
}

#endif // SAMPLE_HH_INCLUDED
