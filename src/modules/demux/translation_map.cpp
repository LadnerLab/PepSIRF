#include "translation_map.h"

namespace codon_aa_mappings
{


    std::unordered_map<std::string,char>
        get_default_map()
    {
        return std::unordered_map<std::string,char>
            {
                { "AAA", 'K' },
                { "AAG", 'K' },
                { "AAC", 'N' },
                { "AAT", 'N' },
                { "ACA", 'T' },
                { "ACC", 'T' },
                { "ACG", 'T' },
                { "ACT", 'T' },
                { "ATA", 'I' },
                { "ATC", 'I' },
                { "ATG", 'M' },
                { "ATT", 'I' },
                { "CAA", 'Q' },
                { "CAC", 'H' },
                { "CAG", 'Q' },
                { "CAT", 'H' },
                { "CCA", 'P' },
                { "CCC", 'P' },
                { "CCG", 'P' },
                { "CCT", 'P' },
                { "CGA", 'R' },
                { "CGC", 'R' },
                { "CGG", 'R' },
                { "CGT", 'R' },
                { "CTA", 'L' },
                { "CTC", 'L' },
                { "CTG", 'L' },
                { "CTT", 'L' },
                { "GAA", 'E' },
                { "GAG", 'E' },
                { "GAC", 'D' },
                { "GAT", 'D' },
                { "GCA", 'A' },
                { "GCG", 'A' },
                { "GCT", 'A' },
                { "GCC", 'A' },
                { "GGA", 'G' },
                { "GGG", 'G' },
                { "GGC", 'G' },
                { "GGT", 'G' },
                { "GTA", 'V' },
                { "GTC", 'V' },
                { "GTG", 'V' },
                { "GTT", 'V' },
                { "TAC", 'Y' },
                { "TAT", 'Y' },
                { "AGT", 'S' },
                { "AGC", 'S' },
                { "AGA", 'R' },
                { "AGG", 'R' },
                { "TCA", 'S' },
                { "TCC", 'S' },
                { "TCG", 'S' },
                { "TCT", 'S' },
                { "TGC", 'C' },
                { "TGG", 'W' },
                { "TGT", 'C' },
                { "TTA", 'L' },
                { "TTC", 'F' },
                { "TTG", 'L' },
                { "TTT", 'F' },
                // begin stop codons
                { "TAG", '_' },
                { "TAA", '_' },
                { "TGA", '_' }
            };
    }

}; // namespace codon_aa_mappings
