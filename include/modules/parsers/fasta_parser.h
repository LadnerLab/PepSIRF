#ifndef FASTA_PARSER_HH_INCLUDED
#define FASTA_PARSER_HH_INCLUDED

#include <stdexcept>
#include <fstream>
#include <iostream>
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
     * @note The '>' character is removed from the name of each sequence, so '>Seq 1' becomes
     *       'Seq 1'.
     * @returns vector of sequences, one for each sequence found in 'filename'
     **/
    std::vector<sequence> parse( std::string filename );

};
#endif // FASTA_PARSER_HH_INCLUDED

