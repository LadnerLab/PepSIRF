#ifndef NT_AA_TRANSLATOR_INCLUDED
#define NT_AA_TRANSLATOR_INCLUDED
#include <string>
#include <assert.h>
#include "translation_map.h"
#include "sequence.h"

/**
 * Class used to transform sequences of nucleotides into
 * their equivalent AA sequences, as defined by a 
 * translation map.
 * @tparam An NT_AA_Map that will translate string codons into 
 *         character amino acids.
 **/
template<typename NT_AA_Map>
class nt_aa_translator
{
public:

    /**
     * The size of a codon, in nucleotides
     **/
    constexpr static std::size_t CODON_SIZE = 3;

    /**
     * Constructor initializer, using a translation map that 
     * defines how codons translate into amino acids.
     * @param translation_map A map whose operator() takes a 
     *        string codon, and returns a character amino acid.
     **/ 
    nt_aa_translator( const NT_AA_Map& translation_map )
        : map{ translation_map } {}
        
    /**
     * Translate a nucleotide sequence into an amino acid sequence.
     * Uses the map object this was initialized with to define the mapping.
     * @param in_seq The NT sequence to translate
     * @returns an AA sequence.
     **/
    sequence translate( const sequence& in_seq ) const
    {
        assert( !( in_seq.seq.size() % CODON_SIZE ) );
        std::string translated;
        translated.reserve( in_seq.seq.size() / CODON_SIZE );

        for( std::string::size_type codon_idx = 0;
             codon_idx < in_seq.seq.size();
             codon_idx += CODON_SIZE
             )
            {
                typename NT_AA_Map
                    ::AAType aa = map( in_seq.seq.substr( codon_idx,
                                                          CODON_SIZE
                                                          )
                                       ); 
                translated.push_back( aa );
            }

        return sequence{ in_seq.name, translated };
    }

    /**
     * Translate a NT sequence into an AA sequence.
     * @note equivalent to calling this->translate( in_seq )
     * @returns an AA sequence.
     **/
    sequence operator()( const sequence& in_seq ) const
    {
        return translate( in_seq );
    }

private:

    /**
     * The map that is used by this object
     **/
    const NT_AA_Map map;

};

#endif 
