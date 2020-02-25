#ifndef NT_AA_TRANSLATOR_INCLUDED
#define NT_AA_TRANSLATOR_INCLUDED
#include <string>
#include <assert.h>
#include "translation_map.h"
#include "sequence.h"

template<typename NT_AA_Map>
class nt_aa_translator
{
public:
    constexpr static std::size_t CODON_SIZE = 3;

    nt_aa_translator( const NT_AA_Map& translation_map )
        : map{ translation_map } {}
        
    sequence operator()( const sequence& in_seq ) const
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

private:
    const NT_AA_Map map;

};

#endif 
