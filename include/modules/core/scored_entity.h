#ifndef SCORED_ENTITY_HH_INCLUDED
#define SCORED_ENTITY_HH_INCLUDED

/**
 * A scored quantity is a  
 **/
template<typename Scored, typename ScoreType>
    class scored_entity
{
 public:
    scored_entity() = default;

 private:
    Scored score_key;
    ScoreType score;


}


#endif // SCORED_ENTITY_HH_INCLUDED
