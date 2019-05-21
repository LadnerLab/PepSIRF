#ifndef FASTA_PARSER_HH_INCLUDED
#define FASTA_PARSER_HH_INCLUDED

#include <stdexcept>
#include <iostream>
#include <fstream>
#include <iterator>
#include <vector>

#include "sequence.h"

/**
 * Parse a fasta file. Given a string filename, parse the sequences in that filename.
 * Parsed data is stored as a vector of 'sequence' classes.
 **/
class fasta_parser
{
 public:
    fasta_parser(); //!< Default constructor.

    /**
     * Parse a fasta file. Records are stored in members of the sequence class.
     * @param filename Name of file to parse.
     * @returns A vector containing sequences, one entry for each '>' character found in the file.
     * @throws Exception upon filename not a valid file.
     **/
    std::vector<sequence> parse( std::string filename );
    


};
#endif // FASTA_PARSER_HH_INCLUDED

