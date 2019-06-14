#ifndef SEQUENCE_INDEXER_HH_INCLUDED
#define SEQUENCE_INDEXER_HH_INCLUDED
#include <vector>
#include <algorithm>

#include "sequence.h"

/**
 * Class for indexing sequences. Takes a set of 
 * reference sequences and sorts them in order of hamming distance 
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
     * but you want the length of the origin sequence to be exactly the length 
     * of all of the sequences in seqs. 
     * @param seqs Reference to vector holding sequences to be copied.
     * @param origin The 'origin' sequence. This origin is analogous to the point 
     *        (0,0) on the xy plane. Each sequence in seqs is ordered in terms of its distance
     *        from this origin sequence. 
     * @note The index is managed internally, to query this index call sequence_indexer::index
     * @note All sequences must be exactly the same length, as hamming distance is not defined for 
     *       strings of different length.
     **/
    void index( std::vector<sequence>& seqs, std::string& origin );

    /**
     * Query the indexed sequences to find all sequences that are 
     * similar to query_seq. For two sequences a and b, and threshold c we say a and b are similar 
     * iff the hamming distance D between a and b is less than or equal to c, i.e.
     * a is similar to b iff D( a, b ) <= c. Pointers to sequences that are 
     * found to be similar to query_seq are placed into the 
     * results vector as a pair, where the second item of the pair is the integer distance 
     * of the sequence to the query_seq.
     * @param results Reference to a vector whose entries are ordered pairs.
     *        The first value in each pair is a pointer to the sequence that was similar 
     *        to query_seq.
     * @param query_seq Reference to sequence that we are searching for the sequences similar to.
     * @param max_dist The maximum allowable distance between sequences in order for them to be 
     *        considered similar.
     * @returns The number of sequences that were found similar to query_seq.
     * @note This method is thread safe with respect to the internal data structures of 
     *       sequence_indexer, but not with respect to the results vector.
     **/
    unsigned int query( std::vector<std::pair<sequence*,int>>& results,
                        sequence& query_seq, int max_dist
                      );

    unsigned int query( std::vector<std::pair<sequence*,int>>& results,
                        sequence& query_seq, int max_dist,
                        bool normalize_sizes
                      );


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

    std::string origin_seq;

 private:
    /**
     * Calcualte the edit distance D between strings s1 and s2. 
     * @param s1 The first string to compare to.
     * @param s2 The second string to compare to.
     * @returns the integer hamming distance between s1 and s2.
     **/
    int edit_distance( const std::string& s1, const std::string& s2 );
};

/**
 * Comparison operator given two nodes.
 **/
bool operator<( sequence_indexer::node const& n1,
                sequence_indexer::node const& n2
              );

#endif // SEQUENCE_INDEXER_HH_INCLUDED

