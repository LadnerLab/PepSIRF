#ifndef PEPTIDE_HH_INCLUDED
#define PEPTIDE_HH_INCLUDED

#include <string>

/**
 * A peptide is a continuous stretch of AA or NT.
 * It can have a name and a sequence.
 **/
class peptide
{
 public:

    /**
     * Default initializer
     **/
    peptide() = default;

    /**
     * Complete initializer for the peptide.
     * @param name The name of the peptide.
     * @param sequence The AA or NT peptide sequence.
     **/
    peptide( const std::string& name,
             const std::string& sequence
           ) : name( name ),  sequence( sequence ) {}

    /**
     * Initialize a peptide with just a sequence.
     * @param sequence The NT or AA sequence to initialize this class with
     **/
    peptide( const std::string& sequence )
        : peptide( "", sequence ) {}


    /**
     * Return a constant reference to the 
     * peptide's name.
     **/
    const std::string& get_name() const;

    /**
     * Return a constant reference to the peptide's peptide.
     **/
    const std::string& get_sequence() const;

    /**
     * Set the peptide's name.
     * @param new_val The new name for this peptide.
     *        This value will be copied into the peptide's 'name'
     *        member. Calling get_name() after set_name() will 
     *        return a string with the contents of new_val
     **/
    void set_name( const std::string& new_val );

    /**
     * Set the peptide's sequence.
     * @param new_val The new sequence for this peptide.
     *        This value will be copied into the peptide's 'sequence'
     *        member. Calling get_sequence() after set_sequence() will 
     *        return a string with the contents of new_val
     **/
    void set_sequence( const std::string& new_val );

    /**
     * Determine whether this peptide equals another.
     * Two peptides are equal if their sequence members 
     * are equal.
     * @param other The peptide to compare with.
     * @returns true if this == other,
     *          false otherwise.
     **/
    bool operator==( const peptide& other ) const;

    /**
     * Determine if this peptide is NOT equal to 
     * another peptide.
     * @param other The peptide to compare with.
     * @note this is equivalent to calling !(a==b)
     * @returns boolean false if this == other, true otherwise
     **/
    bool operator!=( const peptide& other ) const;

 protected:

    /**
     * The peptide's name. 
     * This can be used to uniquely identify a peptide.
     **/
    std::string name;

    /**
     * The peptide's sequence, i.e. the peptide itself.
     **/
    std::string sequence;

};

#endif // PEPTIDE_HH_INCLUDED
