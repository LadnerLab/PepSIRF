#ifndef KMER_TOOLS_HH_INCLUDED
#define KMER_TOOLS_HH_INCLUDED
#include "scored_entity.h"

namespace kmer_tools
{
    constexpr char UNKNOWN_AMINO_ACID = 'X';

    template<template<class, class...> class Dtype>
        int get_kmers( Dtype<std::string>& kmers,
                   std::string seq,
                   int k
                 )
    {
        std::size_t num_kmers = ( seq.length() - k ) + 1;
        kmers.reserve( num_kmers );

        for( std::size_t index = 0; index < num_kmers; ++index )
            {
                std::string kmer = seq.substr( index, k );
                if( kmer.find( 'X' ) == std::string::npos )
                    {
                        kmers.insert( kmers.end(), kmer );
                    }
            }

        return kmers.size();
    }

    template<
    template<typename...> class Set>
        void get_kmer_frequencies( Set<scored_entity<std::string,size_t>>&
                               kmer_frequencies,
                               const std::vector<sequence>& proteins,
                               const std::size_t k
                             )
    {
        std::vector<std::string> current_kmers;

        for( const auto& sequence : proteins )
            {
                get_kmers( current_kmers,
                           sequence.seq,
                           k
                         );

                for( const auto& kmer : current_kmers )
                    {
                        scored_entity<std::string,size_t> score( kmer, 0 );

                        if( kmer_frequencies.find( score )
                            == kmer_frequencies.end()
                          )
                            {
                                kmer_frequencies.emplace( score );
                            }

                        ++kmer_frequencies.find( score )->get_score();
                    }

                current_kmers.clear();
            }

    }

}; // namespace kmer_tools



#endif // KMER_TOOLS_HH_INCLUDED
