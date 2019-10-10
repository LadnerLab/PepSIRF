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

 private:
    std::string name;
    std::string sequence;

};

#endif // PEPTIDE_HH_INCLUDED
