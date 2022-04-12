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
#include "fif_parser.h"

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
    std::vector<sample> parse( const options_demux *d_opts, const std::vector<flex_idx> flexible_idx_data );

    /**
     * Verifies the sample IDs in the samplelist columns used are found in the index id set. Prints missing IDs as warning.
     * @param index_seq_ids the unordered set that is compared against by the sample ids to verify ids used exist.
     * @param sample_ids the unordered set that is checked to verify all ids used in the samplelist are existent in the index/barcode list.
     * 
     * @returns void
     **/
    void check_samples( std::unordered_set<std::string>& index_seq_ids, std::unordered_set<std::string>&  sample_ids );
};

#endif // SAMPLELIST_PARSER_HH_INCLUDED
