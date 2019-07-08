#ifndef FASTQ_SEQUENCE_HH_INCLUDED
#define FASTQ_SEQUENCE_HH_INCLUDED

#include <vector>
#include <string>
#include <assert.h>

#include "sequence.h"

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

class fastq_sequence : public sequence
{

 public:

    fastq_sequence( const std::string& in_name,
                    const std::string& in_seq,
                    const std::string& score_str
                  ) :
                  sequence( in_name, in_seq )
        {

            base = phred::phred_base::PHRED33;
        }

    fastq_sequence( const std::string& in_name,
                    const std::string& in_seq,
                    const std::string& score_str,
                    int ascii_base
                  ) :
                  sequence( in_name, in_seq )
        {

            base = ascii_base;

            // make sure one of the standard bases is used
            assert( base == phred::phred_base::PHRED33
                    || base == phred::phred_base::PHRED64
                  );
        }

    fastq_sequence( const std::string& in_name,
                    const std::string& in_seq
                  ) :
                  sequence( in_name, in_seq )
        {
            base = 0;
        }

    std::vector<unsigned char> score_vec;
    int base;

};



#endif // FASTQ_SEQUENCE_HH_INCLUDED
