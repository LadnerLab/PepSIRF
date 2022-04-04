#ifndef FIF_PARSER_HH_INCLUDED
#define FIF_PARSER_HH_INCLUDED
#include <string>
#include <vector>
#include <unordered_set>
#include <stdexcept>
#include <fstream>
#include <iostream>
#include <utility>
#include <boost/algorithm/string.hpp>

class flex_idx
    {
        public:
            std::string idx_name; // name corresponding to header name in samplesheet
            std::string read_name; // name should either be R1 or R2
            std::size_t idx_start; // index start location (0-based, inclusive)
            std::size_t idx_len; // index length
            std::size_t num_mismatch; // number of mismatches allowed
            std::unordered_set<std::string> barcode_ids;

            flex_idx( std::string col1, std::string col2, std::string col3, std::string col4, std::string col5 )
                {
                    idx_name = col1;
                    read_name = col2;
                    idx_start = std::stoi( col3 );
                    idx_len = std::stoi( col4 );
                    num_mismatch = std::stoi( col5 );
                }

            flex_idx( std::string col1, std::string col2, std::size_t col3, std::size_t col4, std::size_t col5 )
                {
                    idx_name = col1;
                    read_name = col2;
                    idx_start = col3 ;
                    idx_len = col4;
                    num_mismatch = col5;
                }
    };
class fif_parser
    {
        public:
            /**
             * Parse flexible index file containing index information.
             * See help message for (--fif,-f) for more info.
             * @param d_opts contains name of input file to generate flexible index
             *               data and useful info to identify index data for return vec.
             * @return returns vector with each element as data for a single index.
            */
            std::vector<flex_idx> parse( const std::string fif_fname );

    };

#endif // FIF_PARSER_HH_INCLUDED
