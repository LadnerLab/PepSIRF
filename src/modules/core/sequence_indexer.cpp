#include <cstring>

#include "sequence_indexer.h"

sequence_indexer::sequence_indexer() = default;

void sequence_indexer::index( std::vector<sequence>& seqs, std::string& origin )
{
    std::size_t index = 0;
    int distance = 0;

    origin_seq = origin;

    indexed_seqs.reserve( seqs.size() );

    for( index = 0; index < seqs.size(); ++index )
        {
            distance = edit_distance( origin, seqs[ index ].seq );

            indexed_seqs.emplace_back( distance, &seqs[ index ] );
        }
    std::sort( indexed_seqs.begin(), indexed_seqs.end() );
}
unsigned int sequence_indexer::query( std::vector<std::pair<sequence*,int>>& results,
                                      sequence& query_seq, int max_dist,
                                      bool normalize_sizes
                                    )
{
    if( normalize_sizes )
        {
            return query( results, query_seq, max_dist );
        }
    return query( results, query_seq, max_dist + _size_difference( query_seq.seq, origin_seq ) );
}

unsigned int sequence_indexer::query( std::vector<std::pair<sequence*,int>>& results,
                                      sequence& query_seq, int max_dist
                                    )
{
    int orig_distance    = 0;
    int current_distance = 0;
    int lev_distance     = 0;
    unsigned int matches = 0;

    std::vector<node>::iterator it;

    orig_distance = edit_distance( query_seq.seq, origin_seq );

    // skip over the sequences whose distance is greater than
    // our distance + max_dist, we don't need to compute levenshtein distance
    it = std::upper_bound( indexed_seqs.begin(), indexed_seqs.end(),
                           node( orig_distance - max_dist - 1, nullptr )
                         );
    current_distance = it->distance;

    while( std::abs( current_distance - orig_distance ) <= max_dist
           && it != indexed_seqs.end()
         )
        {
            lev_distance = edit_distance( it->seq->seq, query_seq.seq );

            if( lev_distance <= max_dist)
                {
                    results.emplace_back( std::make_pair( it->seq,
                                                          lev_distance
                                                        )
                                        );
                    ++matches;
                }
            ++it;
            current_distance = it->distance;
        }
    return matches;
}


sequence_indexer::node::node( int in_dist, sequence *in_seq )
{
    distance = in_dist;
    seq      = in_seq;
}

sequence_indexer::node::node() = default;

bool sequence_indexer::node::operator<( const node& compare )
{
    return distance < compare.distance;
}

int sequence_indexer::edit_distance( const std::string& s1, const std::string& s2 )
{
    int distance = 0;
    std::size_t index = 0;

    for( index = 0; index < s1.length(); ++index )
        {
            distance += s1[ index ] != s2[ index ];
        }

    return distance;
}

bool operator<( sequence_indexer::node const& n1,
                sequence_indexer::node const& n2
              )
{
    return n1.distance < n2.distance;
}

inline std::size_t
sequence_indexer::_size_difference( const std::string& s1,
                                    const std::string& s2
                                  ) const
{
    std::size_t s1_l = s1.length();
    std::size_t s2_l = s2.length();
    return  s1_l > s2_l ? s1_l - s2_l : s2_l - s1_l;
}
