#ifndef FASTQ_SEQUENCE_HH_INCLUDED
#define FASTQ_SEQUENCE_HH_INCLUDED

#include <vector>
#include <string>

#include "sequence.h"

class fastq_sequence : public sequence
{

 public:

    fastq_sequence( const std::string& in_name,
                    const std::string& in_seq,
                    const std::string& score_str
                  ) :
                  sequence( in_name, in_seq )
        {

        }

    std::vector<unsigned char> score_vec;
};

namespace phred
{

    constexpr int PHRED33_ASCII_BASE = 33;
    constexpr int PHRED64_ASCII_BASE = 64;

    enum phred_base
    {
        PHRED33 = PHRED33_ASCII_BASE,
        PHRED64 = PHRED64_ASCII_BASE
    };

};


#endif // FASTQ_SEQUENCE_HH_INCLUDED
