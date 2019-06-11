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
     *        one entry in the fastq file. So 4 * max_num_records lines will be read from the file, and at most 
     *        max_num_records entries will be added to seq_vector. If max_num_records is set to 0, the entire 
     *        file will be read.
     * @returns bool true if any records were read from the file, false if zero records were 
     *          read from the file.
     **/
    bool parse( std::ifstream& input_file, std::vector<sequence>& seq_vector, std::size_t max_num_records );

    /**
     * Parse an entire fastq file at once, store the parsed sequence data  in 
     * a vector of sequences.
     * @param fname Name of the file to parse.
     * @returns a vector of sequences, one entry in the vector per sequence in the file.
     * @note a call to this function is equivalent to creating an ifstream, a sequence vector and 
     *       calling fastq_parser::parse( ifstream, seq_vector, 0 )
     **/
    std::vector<sequence> parse( const std::string fname );

 private:
    /**
     * Enum member denoting where each item of a fastq entry will 
     * be located if read into an array.
     **/
    enum fastq_item_indices { SEQUENCE_NAME = 0,
                              SEQUENCE,
                              SEQUENCE_REPEAT,
                              SEQUENCE_SCORE
                            };
          

};


#endif // FASTQ_PARSER_HH_INCLUDED
