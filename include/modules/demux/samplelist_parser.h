#ifndef SAMPLELIST_PARSER_HH_INCLUDED
#define SAMPLELIST_PARSER_HH_INCLUDED
#include <string>
#include <vector>
#include <stdexcept>
#include <fstream>
#include <utility>

#include "sample.h" 

class samplelist_parser
{
 public:
    /**
     * Parse a tab-delimited file containing samples, one per line.
     * @param filename The name of the file to parse.
     * @returns vector of samples, one per line in the input file.
     **/
    std::vector<sample> parse( const std::string filename );
};

#endif // SAMPLELIST_PARSER_HH_INCLUDED
