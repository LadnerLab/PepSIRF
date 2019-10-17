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
     **/
    ScoreType score;


};


#endif // SCORED_ENTITY_HH_INCLUDED
