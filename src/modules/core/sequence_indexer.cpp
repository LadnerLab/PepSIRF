#include <cstring>

#include "sequence_indexer.h"
#include "edlib.h"

sequence_indexer::sequence_indexer() = default;

void sequence_indexer::index( std::vector<sequence>& seqs, std::string& origin )
{
    std::size_t index = 0;
    std::size_t size_origin = origin.length();

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
