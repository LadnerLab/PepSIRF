#include <cstring>
#include <assert.h>

#include "sequence_indexer.h"

sequence_indexer::sequence_indexer() = default;

void sequence_indexer::index( std::vector<sequence>& seqs )
{
    std::vector<sequence>::iterator seq_iter = seqs.begin();

    while( seq_iter != seqs.end() )
        {
            tree.insert( &(*seq_iter) );
            ++seq_iter;

            assert( seq_iter->seq.length() == seqs.begin()->seq.length() );

            if(  seq_iter->seq.length() == seqs.begin()->seq.length() )
                {
                    throw std::runtime_error( "The samplelist file is not formatted correctly!" );
                }
        }
}

unsigned int sequence_indexer::query( std::vector<std::pair<sequence*,int>>& results,
                                      sequence& query_seq, std::size_t max_dist
                                    )
{
    unsigned int matches = 0;

    auto result = tree.find_within( &query_seq, max_dist );

    for( auto iter = result.begin(); iter < result.end(); ++iter )
        {
            results.emplace_back( std::make_pair( iter->first, iter->second ) );
            std::push_heap( results.begin(), results.end(), cmp_pairs() );
            ++matches;
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

