#ifndef TEST_UTILS_H_INCLUDED
#define TEST_UTILS_H_INCLUDED
#include <string>
#include <sstream>
#include <fstream>
#include <boost/algorithm/string.hpp>

#include "scored_entity.h"


namespace test_utils
{
    template<template<class K, class...> class Set>
        void parse_kmer_frequency( Set<scored_entity<std::string,std::size_t>>&
                                   dest,
                                   std::ifstream& file_stream
                                 )
        {
            std::string line;
            std::vector<std::string> tokens;

            while( std::getline( file_stream, line ) )
                {
                    boost::split( tokens, line, boost::is_any_of( "\t" ) );
                    std::size_t score = boost::lexical_cast<int>( tokens[ 1 ] );

                    dest.insert( scored_entity<std::string,std::size_t>(  tokens[ 0 ], score ) );

                    tokens.clear();
                }
        }
    
}; // namespace test_utils


#endif // TEST_UTILS_H_INCLUDED
