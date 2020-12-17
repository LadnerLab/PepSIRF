#ifndef SAMPLELIST_PARSER_HH_INCLUDED
#define SAMPLELIST_PARSER_HH_INCLUDED
#include <string>
#include <vector>
#include <unordered_set>
#include <stdexcept>
#include <fstream>
#include <utility>
#include <boost/algorithm/string.hpp>

#include "sample.h" 

class samplelist_parser
{
 public:
    /**
     * Parse a tab-delimited file containing samples, one per line.
     * @param filename The name of the file to parse.
     * @param header_names A set of names defining the columns to be used in sampelist.
     * @returns vector of samples, one per line in the input file.
     **/
    std::vector<sample> parse( const std::string filename, std::unordered_set<std::string> header_names );
};

#endif // SAMPLELIST_PARSER_HH_INCLUDED
