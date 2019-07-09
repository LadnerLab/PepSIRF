#define CATCH_CONFIG_MAIN
#include "catch.hpp"

#include <iostream>
#include <fstream>
#include <omp.h>
#include <algorithm>
#include <vector>
#include <unordered_map>

#include "sequence_indexer.h"
#include "sequence.h"
#include "fastq_sequence.h"
#include "maps.h"
#include "fastq_parser.h"
#include "module_demux.h"
#include "samplelist_parser.h"
#include "sample.h"
#include "fastq_score.h"

static std::string fasta_ex = ">Example sequence 1\nDJFKDLFJSF\nDJFKDJFDJKFJKDF\n\nJDSKFLJFDKSFLJ\n>Example Sequence 2\n\n\nDJFKLSFJDKLSFJKSLFJSKDLFJDKSLFJKLDKDKDK\n\n\n";
TEST_CASE( "Sequence", "[sequence]" )
{
    std::string n = "Sequence name";
    std::string s = "Sequence";
    sequence seq( n, s );

    REQUIRE( seq.name.compare( "Sequence name" ) == 0 );
    REQUIRE( seq.seq.compare( "Sequence" ) == 0 );

    std::string s2 = "ATGCGGGGTC";

    seq.set_seq( s2 );

    REQUIRE( seq.seq.compare( s2 ) == 0 );
    REQUIRE( seq.count( 'A' ) == 1 );
    REQUIRE( seq.count( 'Z' ) == 0 );

    s2 = "AAA";
    REQUIRE( seq.seq.compare( s2 ) != 0 );
}

TEST_CASE( "Parse Fasta", "[fasta_parser]" )
{
    fasta_parser fp;
    std::vector<sequence> vec;
    vec = fp.parse( "../test/test.fasta" );

    REQUIRE( vec.size() == 110 );

    unsigned int index = 0;
    for( index = 0; index < vec.size(); ++index )
        {
            REQUIRE( vec[ index ].seq.compare( "" ) );
            REQUIRE( vec[ index ].name.compare( "" ) );
        }

    std::ofstream stream;
    stream.open( "out.fasta" );
    stream << fasta_ex;
    stream.close();

    vec = fp.parse( "out.fasta" );
    REQUIRE( vec.size() == 2 );
    REQUIRE( !vec[ 0 ].name.compare( "Example sequence 1" ) );
    REQUIRE( !vec[ 0 ].seq.compare( "DJFKDLFJSFDJFKDJFDJKFJKDFJDSKFLJFDKSFLJ" ) );
    REQUIRE( !vec[ 1 ].name.compare( "Example Sequence 2" ) );
    REQUIRE( !vec[ 1 ].seq.compare( "DJFKLSFJDKLSFJKSLFJSKDLFJDKSLFJKLDKDKDK" ) );
    REQUIRE( vec[ 1 ].seq.compare( vec[ 0 ].seq ) );
    REQUIRE( vec[ 0 ].seq.compare( vec[ 1 ].seq ) );


    remove( "out.fasta" );

}

TEST_CASE( "Parse Fastq", "[fastq_parser]" )
{
    std::vector<fastq_sequence> seq_vec;
    std::string fastq_fname = "../test/test_fastq.fastq";
    std::ifstream in_file( fastq_fname, std::ios_base::in );

    fastq_parser parse;

    size_t step = 10;

    REQUIRE( parse.parse( in_file, seq_vec, 0 ) == true );

    REQUIRE( seq_vec.size() == 100 );
    size_t index = 0;

    for( index = 0; index < seq_vec.size(); ++index )
        {
            REQUIRE( seq_vec[ index ].seq.compare( "" ) );
            REQUIRE( seq_vec[ index ].name.compare( "" ) );
        }

    // reset file pointer
    in_file.clear();
    in_file.seekg( 0 );

    std::vector<fastq_sequence> seq_vec2;
    seq_vec2.reserve( 100 );

    REQUIRE( seq_vec2.size() == 0 );

    for( index = 0; index < 10; ++index )
        {
            parse.parse( in_file, seq_vec2, step );
        }
    REQUIRE( seq_vec2.size() == 100 );


    std::vector<fastq_sequence> seq_vec3;
    seq_vec3.reserve( 100 );
    // reset file pointer
    in_file.clear();
    in_file.seekg( 0 );

    while( parse.parse( in_file, seq_vec3, 3 ) ) ;

    REQUIRE( seq_vec3.size() == 100 );

    for( index = 0; index < seq_vec3.size(); ++index )
        {
            REQUIRE( seq_vec3[ index ].seq.compare( "" ) );
            REQUIRE( seq_vec3[ index ].name.compare( "" ) );
        }

    for( index = 0; index < seq_vec2.size(); ++index )
        {
            REQUIRE( seq_vec2[ index ].seq.compare( "" ) );
            REQUIRE( seq_vec2[ index ].name.compare( "" ) );
        }
    for( index = 0; index < seq_vec2.size(); ++index )
        {
            REQUIRE( !seq_vec2[ index ].seq.compare( seq_vec[ index ].seq ) );
            REQUIRE( !seq_vec2[ index ].name.compare( seq_vec[ index ].name ) );
            REQUIRE( !seq_vec3[ index ].name.compare( seq_vec[ index ].name ) );
            REQUIRE( !seq_vec3[ index ].seq.compare( seq_vec[ index ].seq ) );
        }

}

TEST_CASE( "Add seqs to map", "[module_demux]" )
{
    fasta_parser fp;
    module_demux mod;

    std::vector<sequence> vec;
    std::size_t num_samples = 12;
    vec = fp.parse( "../test/test.fasta" );

    parallel_map<sequence, std::vector<std::size_t>*> my_map;

    mod.add_seqs_to_map( my_map, vec, num_samples );

    REQUIRE( my_map.size() == 110 );

    parallel_map<sequence, std::vector<std::size_t>*>::iterator it = my_map.begin();
    size_t index = 0;
                       
    while( it != my_map.end() )
        {
            REQUIRE( it->second->size() == num_samples );

            for( index = 0; index < it->second->size(); ++index )
                {
                    REQUIRE( it->second->at( index ) == 0 );
                }
            ++it;
        }

}

TEST_CASE( "Test samplelist_parser and sample", "[samplelist_parser]" )
{
    samplelist_parser slp;
    std::string filename = "../test/test_samplelist.tsv";

    auto vec = slp.parse( filename );

    REQUIRE( vec.size() == 96 );


    unsigned int index = 0;
    sequential_map<sample, int> s_map( vec.size() );
    for( index = 0 ; index < vec.size(); ++index )
        {
            s_map[ vec[ index ] ] = vec[ index ].id;
        }

    REQUIRE( s_map.size() <= vec.size() );
}

TEST_CASE( "Test String Indexing", "[string_indexer]" )
{
    std::string origin = "";
    sequence_indexer si;
    fastq_parser fp;
    auto seqs = fp.parse( "../test/test_fastq.fastq" );
    std::vector<std::pair<sequence*,int>> result_set;

    std::vector<sequence> seqs_fastq;

    std::for_each( seqs.begin(), seqs.end(),
                   [&]( const sequence& seq )
                   {
                       seqs_fastq.emplace_back( seq.name, seq.seq );
                   }
                 );


    si.index( seqs_fastq );
    REQUIRE( si.tree.size() > 0 );
    REQUIRE( si.tree.size() == seqs.size() );

    std::string q1 = "NGCCAGCTTGCGGCAAAACTGCGTAACCGTCTTCTCGTTCTCTAAAAACCATTTTTCGTCCCCTTCGGGGCGGTGGTCTATAGTGTTATTAATATCAAGTTGGGGGAGCACATTGTAGCATTGTGCCAATTCATCCATTAACTTCTCAGT";
    std::string q2 = "CTCCAGCTTGCGGCAAAACTGCGTAACCGTCTTCTCGTTCTCTAAAAACCATTTTTCGTCCCCTTCGGGGCGGTGGTCTATAGTGTTATTAATATCAAGTTGGGGGAGCACATTGTAGCATTGTGCCAATTCATCCATTAACTTATCAGT";
    sequence s1( "", q1 );
    sequence s2( "", q2 );

    unsigned int num_matches = si.query( result_set, s1, 0 );

    REQUIRE( num_matches == 1 );
    REQUIRE( result_set.size() == 1 );
    REQUIRE( std::get<1>(result_set[ 0 ]) == 0 );

    result_set.clear();
    num_matches = si.query( result_set, s2, 3 );

    REQUIRE( num_matches == 1 );
    REQUIRE( result_set.size() == 1 );
    REQUIRE( std::get<1>(result_set[ 0 ]) == 3 );

    result_set.clear();
    num_matches = si.query( result_set, s2, 150 );

    REQUIRE( num_matches == 100 );
    REQUIRE( result_set.size() == 100 );
    REQUIRE( std::get<1>(result_set[ 0 ]) == 3 );

    std::vector<sequence> seqs2;
    seqs2.emplace_back( "", "ATGC" );
    seqs2.emplace_back( "", "TGC" );

    sequence_indexer si2;

}

TEST_CASE( "Test Count Generation", "[module_demux]" )
{

    fastq_parser fp;
    module_demux mod;
    sequence_indexer lib_idx;

    std::vector<fastq_sequence> vec_a;
    std::vector<sequence> vec;
    std::size_t num_samples = 1;

    std::size_t seq_length     = 150;
    std::size_t seq_start      = 0;
    std::size_t num_mismatches = 0;
    
    vec_a = fp.parse( "../test/test_fastq.fastq" );

    std::for_each( vec_a.begin(), vec_a.end(),
                   [&]( const sequence& seq )
                   {
                       vec.emplace_back( seq.name, seq.seq );
                   }
                 );


    lib_idx.index( vec );

    parallel_map<sequence, std::vector<std::size_t>*> my_map;

    mod.add_seqs_to_map( my_map, vec, num_samples );

    parallel_map<sequence, std::vector<std::size_t>*>::iterator it = my_map.begin();
    parallel_map<sequence, std::vector<std::size_t>*>::iterator seq_match;

    std::size_t index = 0;

    // perfect matches
    for( index = 0; index < my_map.size(); ++index )
        {
                seq_match = mod._find_with_shifted_mismatch( my_map, vec[ index ],
                                                             lib_idx, num_mismatches,
                                                             seq_start, seq_length
                                                           );
            REQUIRE( seq_match != my_map.end() );
            REQUIRE( seq_match->second->size() == 1 );

            seq_match->second->at( 0 )++;

            REQUIRE( seq_match->second->at( 0 ) == 1 );

        }
    sequence seq =
    sequence( "", "AGCCAGCTTGCGGCAAAACTGCGTAACCGTCTTCTCGTTCTCTAAAAACCATTTTTCGTCCCCTTCGGGGCGGTGGTCTATAGTGTTATTAATATCAAGTTGGGGGAGCACATTGTAGCATTGTGCCAATTCATCCATTAACTTCTCAGT" );

    seq_match = mod._find_with_shifted_mismatch( my_map, seq,
                                                 lib_idx, num_mismatches,
                                                 seq_start, seq_length
                                                 );
    REQUIRE( seq_match == my_map.end() );

    seq_match = mod._find_with_shifted_mismatch( my_map, seq,
                                                 lib_idx, 1,
                                                 seq_start, seq_length
                                                 );
    REQUIRE( seq_match != my_map.end() );
    seq =
    sequence( "", "AGCCAGCTTGCGGCAAAACTGCGTAACCGTCTTCTCGTTCTCTAAAAACCATTTTTCGTCCCCTTCGGGGCGGTGGTCTATAGTGTTATTAATATCAAGTTGGGGGAGCACATTGTAGCATTGTGCCAATTCATCCATTAACTTCGCAAT" );

    seq_match = mod._find_with_shifted_mismatch( my_map, seq,
                                                 lib_idx, 3,
                                                 seq_start, seq_length
                                                 );
    REQUIRE( seq_match != my_map.end() );




}

TEST_CASE( "Fastq_scorer", "[fastq_score]" )
{

    fastq_parser fp;
    auto vec = fp.parse( "../test/test_fastq.fastq" );

    for( auto it = vec.begin(); it != vec.end(); ++it )
        {
            auto score = fastq_score::get_avg_score( it->scores.begin(),
                                                     it->scores.end(),
                                                     fastq_score::phred::PHRED33
                                                     );
            REQUIRE( score >= 0 );
        }
}
