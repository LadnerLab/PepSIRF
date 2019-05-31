#ifndef FASTQ_PARSER_HH_INCLUDED
#define FASTQ_PARSER_HH_INCLUDED

#include <vector>
#include <fstream>
#include <string>
#include <cstddef>

#include "sequence.h"

class fastq_parser
{
 public:
    /**
     * Default constructor
     **/
    fastq_parser();

    /**
     * Parse at most max_num_records from the ifstream input_file, storing sequences in seq_vector.
     * @param input_file Reference to an open ifstream.
     * @param seq_vector Vector containing (or an empty vector) sequences.
     * @param max_num_records The maximum number of records to parse. Note that a 'record' is considered 
     *        one entry in the fastq entry. So 4 * max_num_records will be read from the file, and at most 
     *        max_num_records entries will be added to seq_vector. If max_num_records is set to 0, the entire 
     *        file will be read.
     * @returns bool true if max_num_records were read from the file, false if less than max_num_records were 
     *          read from the file.
     **/
    bool parse( std::ifstream& input_file, std::vector<sequence>& seq_vector, std::size_t max_num_records );

 private:
    enum fastq_item_indices { SEQUENCE_NAME = 0,
                              SEQUENCE,
                              SEQUENCE_REPEAT,
                              SEQUENCE_SCORE
                            };
          

};


#endif // FASTQ_PARSER_HH_INCLUDED
