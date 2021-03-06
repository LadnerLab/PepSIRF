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

    using score_type = double;
    using probe_type = std::string;
    using rank_track_type = std::unordered_map<score_type,
                            std::vector<probe_type>>;

    /**
     * Default constructor, sets rounding factor to zero.
     **/
    probe_rank()
        : rounding_factor{ 0 } {}

    /**
     * Constructor that sets the rounding factor of 
     * this object. 
     * @param rounding_factor the rounding factor to use when 
     *        ranking probes based upon their scores.
     *        Each score is rounded to the nearest 1/10^x for 
     *        rounding factor x.
     **/
    probe_rank( const std::size_t rounding_factor )
        : rounding_factor{ rounding_factor } {}

    /**
     * Rank a probe given its score.
     * Scores with the same rank (score rounded to the nearest 
     * 1/10^x for rounding_factor x) are stored together.
     * @param score The score of the probe to rank
     * @param probe The probe to rank
     **/
    void rank_probe( const score_type score,
                     const probe_type probe
                   );

    /**
     * Get a mutable iterator to the probes with rank.
     * @note Returns rank_track_type::end if rank is valid
     * @note rounds rank to this::rounding_factor before looking 
     * @param rank The rank to find the probes of
     * @returns rank_track_type::iterator pointing to the result of 
     *          looking for the rank. If the rank is not found, the result is 
     *          equal to this::get_probe_ranks::end.
     **/
    rank_track_type::iterator
        get_probes_with_rank( const score_type rank );

    /**
     * Get a const iterator to the probes with rank.
     * @note Returns rank_track_type::end if rank is valid
     * @note rounds rank to this::rounding_factor before looking 
     * @param rank The rank to find the probes of
     * @returns rank_track_type::const_iterator pointing to the result of 
     *          looking for the rank. If the rank is not found, the result is 
     *          equal to this::get_probe_ranks::end.
     **/
    rank_track_type::const_iterator
        get_probes_with_rank( const score_type rank ) const;


    /**
     * Get a constant reference to this object's probe ranks.
     **/
    const rank_track_type& get_probe_ranks() const;

 private:
    rank_track_type ranked_probes;
    const std::size_t rounding_factor;

    double round_to_factor( const double value ) const;
};

#endif // PROBE_RANK_HH_INCLUDED
