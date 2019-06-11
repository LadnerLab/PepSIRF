#include <cstring>

#include "sequence_indexer.h"
#include "edlib.h"

sequence_indexer::sequence_indexer() = default;

void sequence_indexer::index( std::vector<sequence>& seqs, std::string& origin )
{
    std::size_t index = 0;
    std::size_t size_origin = origin.length();
    origin_seq = origin;

    EdlibAlignConfig config = edlibDefaultAlignConfig();
    EdlibAlignResult result;

    indexed_seqs.reserve( seqs.size() );

    for( index = 0; index < seqs.size(); ++index )
        {
            result = edlibAlign( seqs[ index ].seq.c_str(), seqs[ index ].seq.length(),
                                 origin.c_str(), size_origin, config
                               );

            indexed_seqs.emplace_back( result.editDistance, &seqs[ index ] );

            edlibFreeAlignResult( result );
        }
    std::sort( indexed_seqs.begin(), indexed_seqs.end() );
}

unsigned int sequence_indexer::query( std::vector<std::pair<sequence*,int>>& results,
                                      sequence& query_seq, int max_dist
                                    )
{
    EdlibAlignConfig config = edlibDefaultAlignConfig();
    EdlibAlignResult result;

    int orig_distance    = 0;
    int current_distance = 0;
    int lev_distance     = 0;
    unsigned int matches = 0;

    std::vector<node>::iterator it = indexed_seqs.begin();

    result = edlibAlign( query_seq.seq.c_str(), query_seq.seq.length(),
                         origin_seq.c_str(), origin_seq.length(), config
                       );

    orig_distance = result.editDistance;

    current_distance = it->distance;

    edlibFreeAlignResult( result );

    // skip over the sequences whose distance is greater than
    // our distance + max_dist, we don't need to compute levenshtein distance
    while( std::abs( current_distance - orig_distance ) > max_dist
           && it != indexed_seqs.end()
         )
        {
            ++it;
            current_distance = it->distance;
        }

    while( std::abs( current_distance - orig_distance ) <= max_dist
           && it != indexed_seqs.end()
         )
        {
            result = edlibAlign( it->seq->seq.c_str(), it->seq->seq.length(),
                                 query_seq.seq.c_str(), query_seq.seq.length(),
                                 config
                               );
            lev_distance = result.editDistance;

            edlibFreeAlignResult( result );

            if( lev_distance <= max_dist)
                {
                    results.emplace_back( std::make_pair( it->seq,
                                                          lev_distance
                                                        )
                                        );
                    ++matches;
                }
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
