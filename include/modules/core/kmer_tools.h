#ifndef KMER_TOOLS_HH_INCLUDED
#define KMER_TOOLS_HH_INCLUDED

namespace kmer_tools
{

    int get_kmers( std::vector<std::string>& kmers,
                   std::string seq,
                   int k
                 )
    {
        std::size_t num_kmers = ( seq.length() - k ) + 1;
        kmers.reserve( num_kmers );

        for( std::size_t index = 0; index < num_kmers; ++index )
            {
                kmers.emplace_back(
                                   seq.substr( index, k )
                                  );
            }

        return kmers.size();
    }

}; // namespace kmer_tools



#endif // KMER_TOOLS_HH_INCLUDED
