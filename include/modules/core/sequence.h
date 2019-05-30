#ifndef SEQUENCE_HH_INCLUDED
#define SEQUENCE_HH_INCLUDED
#include <string>

class sequence
{

 public:

    /**
     * Default constructor, initializes the class's name and sequence.
     **/
    sequence();

    /**
     * Constructor where the name and sequence are passed in. 
     * @param in_name reference to string name. The sequence's name will 
     *        be set to the value of this string.
     * @param in_seq The genetic data of a sequence. The sequence's name 
     *        will be set to the value of this string.
     **/
    sequence( std::string& in_name, std::string& in_seq );

    /** 
     * Default destructor.
     **/
    ~sequence();

    /**
     * The name of the sequence represented as a std::string.
     **/
    std::string name;

    /**
     * The sequence itself, represented by a string.
     **/
    std::string seq;

    /**
     * Set the name of the sequence to the specified string value.
     * @param new_name Reference to a string, this object's new name will 
     *        be the value of new_name.
     **/
    void set_name( std::string& new_name );

    /**
     * Set the sequence data of the sequence.
     * @param new_seq reference to a string, this object's  new 
     *        name will be the value of new_seq.
     **/
    void set_seq( std::string& new_seq );


    /**
     * Count occurrences of the character val in the sequence.
     * @param val the character to search for in the sequence
     * @returns integer count of the number of times 'val' appears in the sequence.
     **/
    int count( const char val );

};
#endif //SEQUENCE_HH_INCLUDED 
