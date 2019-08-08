#ifndef KMER_TOOLS_HH_INCLUDED
#define KMER_TOOLS_HH_INCLUDED

namespace kmer_tools
{

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
                kmers.insert( kmers.end(), seq.substr( index, k ) );
            }

        return kmers.size();
    }

}; // namespace kmer_tools



#endif // KMER_TOOLS_HH_INCLUDED
