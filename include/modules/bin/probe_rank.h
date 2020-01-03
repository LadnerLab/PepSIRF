#ifndef PROBE_RANK_HH_INCLUDED
#define PROBE_RANK_HH_INCLUDED
#include <string>
#include <unordered_map>
#include <vector>

/**
 * A class to rank probes given their scores. 
 **/
class probe_rank
{
 public:

    /**
     * Default constructor, sets rounding factor to zero.
     **/
    probe_rank()
        : rounding_factor{ 0 } {}

    /**
     * Constructor that sets the rounding factor of 
     * this object. 
     **/
    probe_rank( const std::size_t rounding_factor )
        : rounding_factor{ rounding_factor } {}

    using key_type = double;
    using val_type = std::string;
    using rank_track_type = std::unordered_map<key_type,
                            std::vector<val_type>>;

 private:
    rank_track_type ranked_probes;
    std::size_t rounding_factor;

    double round_to_factor( const double value ) const;
};

#endif // PROBE_RANK_HH_INCLUDED
