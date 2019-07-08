#ifndef FASTQ_SEQUENCE_HH_INCLUDED
#define FASTQ_SEQUENCE_HH_INCLUDED

#include <vector>
#include <string>
#include <assert.h>

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

            scores = score_str;
        }

    fastq_sequence( const std::string& in_name,
                    const std::string& in_seq,
                    const std::string& score_str,
                    int ascii_base
                  ) :
                  sequence( in_name, in_seq )
        {

            base = ascii_base;
            scores = score_str;
        }

    fastq_sequence( const std::string& in_name,
                    const std::string& in_seq
                  ) :
                  sequence( in_name, in_seq )
        {
            base = 0;
        }

    std::string scores;
    int base;

};



#endif // FASTQ_SEQUENCE_HH_INCLUDED
