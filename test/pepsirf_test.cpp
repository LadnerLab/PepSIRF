#define CATCH_CONFIG_MAIN
#include "catch.hpp"

#include <iostream>
#include <fstream>
#include <vector>
#include "sequence.h"
#include "fasta_parser.h"
#include "fastq_parser.h"

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
    REQUIRE( !vec[ 0 ].name.compare( ">Example sequence 1" ) );
    REQUIRE( !vec[ 0 ].seq.compare( "DJFKDLFJSFDJFKDJFDJKFJKDFJDSKFLJFDKSFLJ" ) );
    REQUIRE( !vec[ 1 ].name.compare( ">Example Sequence 2" ) );
    REQUIRE( !vec[ 1 ].seq.compare( "DJFKLSFJDKLSFJKSLFJSKDLFJDKSLFJKLDKDKDK" ) );
    REQUIRE( vec[ 1 ].seq.compare( vec[ 0 ].seq ) );
    REQUIRE( vec[ 0 ].seq.compare( vec[ 1 ].seq ) );


    remove( "out.fasta" );

}

TEST_CASE( "Parse Fastq", "[fastq_parser]" )
{
    std::vector<sequence> seq_vec;
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

    std::vector<sequence> seq_vec2;
    seq_vec2.reserve( 100 );

    REQUIRE( seq_vec2.size() == 0 );

    for( index = 0; index < 10; ++index )
        {
            parse.parse( in_file, seq_vec2, step );
        }
    REQUIRE( seq_vec2.size() == 100 );


    std::vector<sequence> seq_vec3;
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
