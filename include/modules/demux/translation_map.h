#ifndef TRANSLATION_MAP_HH_INCLUDED
#define TRANSLATION_MAP_HH_INCLUDED
#include <unordered_map>
#include <string>

/**
 * A class that is used to associate codons with the amino acids they translate 
 * to. 
 * @tparam Map A map fulfilling the interface of std::map that allows 
 *         us to call map[ codon ]. The key_type of this map should be 
 *         the type that represents a codon, and the mapped_type should be 
 *         a type representing an Amino Acid.
 **/
template<typename Map>
class codon_aa_map
{
 public:

    /**
     * The map type that is used by this translator.
     **/
    using MapType = Map;

    /**
     * A type representing a codon.
     **/
    using CodonType = typename MapType::key_type;

    /**
     * Represents an amino acid.
     **/
    using AAType    = typename MapType::mapped_type;

    /**
     * Construct the codon->aa map.
     * @param map The map defining which codons
     *        translate to which amino acids.
     * @note makes a copy of the map passed as an argument.
     **/
    codon_aa_map( const MapType& map )
        : translator{ map }
    {}

    /**
     * Translate a codon into an amino acid, using the 
     * map stored by this translation map.
     * @param codon the codon to translate into an amino acid.
     * @returns the amino acid corresponding to the specified codon.
     **/
    AAType translate_codon( const CodonType codon ) const
    {
        return translator.at( codon );
    }

    /**
     * Translate a codon into its equivalent amino acid.
     * @note Equivalent to calling this->translate_codon( codon )
     * @returns the amino acid corresponding to the specified codon.
     **/
    AAType operator()( const CodonType codon ) const
    {
        return translate_codon( codon );
    }


 private:

    /**
     * The map that is used to translate codons into amino acids.
     **/
    const MapType translator;
    
};

namespace codon_aa_mappings
{
    /**
     * Get the default Codon -> AA map.
     * This simply constructs and returns the map.
     * This is necessary for the reasons listed in 
     * this post:
     * https://stackoverflow.com/questions/9092479/why-isnt-my-extern-variable-initialized-yet
     **/
    std::unordered_map<std::string,char>
        get_default_map();
    
    /**
     * The default codon_aa map.
     **/
    static const codon_aa_map<decltype( get_default_map() )>
        default_codon_aa_map{ get_default_map() };

    /**
     * The type of the default codon aa map.
     **/
    using default_map_type = decltype( default_codon_aa_map );

}; // namespace codon_aa_mappings

#endif // TRANSLATION_MAP_HH_INCLUDED
