#ifndef SEQUENCE_INDEXER_HH_INCLUDED
#define SEQUENCE_INDEXER_HH_INCLUDED
#include <vector>
#include<algorithm>

#include "sequence.h"

/**
 * Class for indexing sequences. Takes a set of 
 * reference sequences and sorts them in order of Levenshtein distance 
 * to an 'origin' sequence. Once indexed, the sequence_indexer can be 
 * queried (sequence_indexer::query) to find all of the sequences that 
 * are within a maximum distance of the query sequence. By indexing we 
 * are able to reduce the number of comparisons from O(n) to roughly 
 * O(logn) for a single query.
 **/
class sequence_indexer
{
 public:
    /**
     * Default constructor.
     **/
    sequence_indexer();

    /**
     * Index a vector of sequences. Sequences are placed in order of their 
     * distance to the origin sequence. Generally the origin sequence does not matter, 
     * but you want the length of the origin sequence to be close (or exactly) to the length 
     * of all of the sequences in seqs. 
     * @param seqs Reference to vector holding sequences to be copied.
     * @param origin The 'origin' sequence. This origin is analogous to the point 
     *        (0,0) on the xy plane. Each sequence in seqs is ordered in terms of its distance
     *        from this origin sequence. 
     * @note The index is managed internally, to query this index call sequence_indexer::index
     **/
    void index( std::vector<sequence>& seqs, std::string& origin );

    /**
     * A node class, these nodes are stored in a vector used internally by 
     * the sequence_indexer.
     **/
    class node
    {
    public:
        /**
         * Argument constructor. Called by 'emplace_back'.
         * @param in_dist Distance of this sequence from the origin.
         * @param in_seq pointer to a sequence that this node considers.
         **/
        node( int in_dist, sequence *in_seq );

        /**
         * Default constructor.
         **/
        node();

        /**
         * Comparison less than operator. For two nodes a and b, we say 
         * a < b iff a.distance < b.distance.
         * @param compare Node to compare against
         **/
        bool operator<( const node& compare );

        /**
         * Integer distance from the origin. 
         * @note this is a signed integer because the value returned by 
         *       edlib may be < 0.
         **/
        int distance;

        /**
         * Sequence referenced by this node.
         **/
        sequence *seq;
    };

    /**
     * Vector containing nodes. After a call to sequence_indexer::index,
     * this vector will contain a set of nodes sorted on distance from the 
     * origin.
     **/
    std::vector<node> indexed_seqs;
};



#endif // SEQUENCE_INDEXER_HH_INCLUDED
