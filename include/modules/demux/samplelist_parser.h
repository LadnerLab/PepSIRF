#ifndef SAMPLELIST_PARSER_HH_INCLUDED
#define SAMPLELIST_PARSER_HH_INCLUDED
#include <string>
#include <vector>
#include <unordered_set>
#include <stdexcept>
#include <fstream>
#include <iostream>
#include <utility>
#include <boost/algorithm/string.hpp>
#include <module_demux.h>
#include "sample.h" 

class samplelist_parser
{
 public:
    /**
     * Parse a tab-delimited file containing samples, one per line.
     * @param d_opts The pointer references the demux options. Specifically
     *               the samplelist filename and headers: samplename, 
     *               index1 and index2 are accessed to create samples.
     * @returns vector of samples, one per line in the input file.
     **/
    std::vector<sample> parse( const options_demux *d_opts );
};

#endif // SAMPLELIST_PARSER_HH_INCLUDED
