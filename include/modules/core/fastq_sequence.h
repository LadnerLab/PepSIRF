#ifndef FASTQ_SEQUENCE_HH_INCLUDED
#define FASTQ_SEQUENCE_HH_INCLUDED

#include <vector>
#include <string>
#include <assert.h>

#include "sequence.h"


/**
 * Class to represent a fastq sequence. 
 * The only difference between this and the base sequence class 
 * is the presence of a score string. This string can be used to 
 * determine the quality of reads.
 **/
class fastq_sequence : public sequence
{

 public:

    /**
     * Initializer that includes the score threshold. 
     * @param in_name The name of the sequence that 
     *        has been read from the fastq.
     * @param in_seq The sequence that has been read from the fastq
     * @param score_str The score vector of quality for the read.
     **/
    fastq_sequence( const std::string& in_name,
                    const std::string& in_seq,
                    const std::string& score_str
                  ) :
                  sequence( in_name, in_seq )
        {
            scores = score_str;
        }

    /**
     * Initializer that does NOT the score threshold. 
     * @param in_name The name of the sequence that 
     *        has been read from the fastq.
     * @param in_seq The sequence that has been read from the fastq
     **/
    fastq_sequence( const std::string& in_name,
                    const std::string& in_seq
                  ) :
                  sequence( in_name, in_seq ) {}


    /**
     * Scores vector.
     **/
    std::string scores;

};



#endif // FASTQ_SEQUENCE_HH_INCLUDED
