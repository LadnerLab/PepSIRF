#ifndef SCORED_ENTITY_HH_INCLUDED
#define SCORED_ENTITY_HH_INCLUDED

/**
 * A scored entity is an entity that is associated 
 * with a score. 
 **/
template<typename Scored, typename ScoreType>
    class scored_entity
{
 public:
    /**
     * Default constructor.
     **/
    scored_entity() {};

    /**
     * Argument constructor, build a 
     * scored entity out of a key and its score. 
     * @param key The key to store in this object,
     *        it is the entity that is scored.
     * @score score The score to assign to the 
     *        key. 
     **/
    scored_entity( const Scored& key,
                   const ScoreType& score
                   ) : score_key( key ), score( score ) {}

    /**
     * Set the score of this object to a new value.
     * @param score The new score of this object.
     **/ 
    void set_score( const ScoreType& score )
    {
        score = score;
    }

    /**
     * Set the key of this object.
     * @param new_key The new key of this object, 
     *        which is used to differentiate objects of 
     *        scored_entity type.
     **/
    void set_key( const Scored& new_key )
    {
        score_key = new_key;
    }

    /**
     * Return a constant reference to this item's
     * score key.
     * @returns See above.
     **/
    const Scored& get_key() const
    {
        return score_key;
    }

    /**
     * Return a constant reference to this item's
     * score.
     * @returns See above.
     **/
    ScoreType& get_score() const
    {
        return score;
    }


    /**
     * Determine whether this is equal to other.
     * For scored_entities a and b, we say a == b 
     * iff a.score_key == b.score_key 
     * @param other The scored_entity to compare with.
     * @returns boolean true if this == other, false otherwise
     **/
    bool operator==( const scored_entity& other ) const
    {
        return score_key == other.score_key;
    }

    /**
     * Determine whether this is NOT equal to other.
     * For scored_entities a and b, we say a != b
     * iff !( a == b )
     **/
    bool operator!=( const scored_entity& other ) const
    {
        return !( *this == other );
    }

 private:
    /**
     * The 'score_key' for the entity, this is used to 
     * uniquely identify a scored entity when hashing, 
     * and comparing equality.
     **/
    Scored score_key;

    /**
     * The score of the entity associated with 
     * score_key.
     * The score is marked mutable because changing the score
     * does NOT change the object itself.
     **/
    mutable ScoreType score;


};


namespace std
{
    /**
     * Hash for a scored entity. Here, we define the 
     * hash of a scored as the hash of its key.
     **/
    template<typename K, typename V>
        struct hash<scored_entity<K,V>>
        {
            std::size_t operator()( const scored_entity<K,V>& val ) const
                {
                    return std::hash<K>()( val.get_key() );
                }

        };
}; // namespace std

#endif // SCORED_ENTITY_HH_INCLUDED
