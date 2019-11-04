#define CATCH_CONFIG_MAIN
#include "catch.hpp"

#include <iostream>
#include <stdexcept>
#include <fstream>
#include <omp.h>
#include <algorithm>
#include <vector>
#include <unordered_map>
#include <cmath>
#include <unordered_set>
#include <string>
#include <fstream>
#include <streambuf>
#include <boost/lexical_cast.hpp>

#include "test_utils.h"
#include "overlap_data.h"
#include "distance_tools.h"
#include "sequence_indexer.h"
#include "sequence.h"
#include "fastq_sequence.h"
#include "maps.h"
#include "distance_matrix.h"
#include "fastq_parser.h"
#include "module_demux.h"
#include "module_deconv.h"
#include "samplelist_parser.h"
#include "sample.h"
#include "fastq_score.h"
#include "kmer_tools.h"
#include "fs_tools.h"
#include "peptide.h"
#include "scored_peptide.h"
#include "species_id.h"
#include "scored_entity.h"

using namespace util;

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
    std::unordered_map<sample, int> s_map( vec.size() );
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

TEST_CASE( "get_kmers", "[kmer_scores]" )
{
    std::string sequence = "ABCDEFGHIJKLMONPQRSTUVWYZ";
    std::vector<std::string> kmers;
    int num_kmers = 0;

    num_kmers = kmer_tools::get_kmers( kmers, sequence, 1 );
    REQUIRE( num_kmers == 25 );
    REQUIRE( kmers.size() == 25 );

    kmers.clear();

    for( std::size_t index = 0; index < sequence.length(); ++index  )
        {
            num_kmers = kmer_tools::get_kmers( kmers, sequence, index + 1 );
            REQUIRE( num_kmers == sequence.length() - ( index + 1 ) + 1 );
            REQUIRE( kmers.size() == sequence.length() - ( index + 1 ) + 1 );
            kmers.clear();
        }
                                                    
                                                  
}


TEST_CASE( "get_tie_candidates_integer", "[module_deconv]" )
{
    auto mod = module_deconv();

    std::vector<std::pair<std::string, double>> candidates;
    std::vector<std::pair<std::string, double>> scores;
    double threshold      = 4.0;
    double ovlp_threshold = 4.0;

    scores.emplace_back( "100", 245.0 );
    scores.emplace_back( "101", 105.0 );
    scores.emplace_back( "102", 5.0 );
    scores.emplace_back( "103", 245.0 );
    scores.emplace_back( "104", 135.0 );
    scores.emplace_back( "105", 249.0 );
    std::sort( scores.begin(),
               scores.end(),
               util::compare_pair_non_increasing
               <std::string, double>()
            );

    auto t_type = mod.get_tie_candidates( candidates,
                                          scores,
                                          threshold,
                                          ovlp_threshold,
                                          difference<double>()
                                        );
                           
    REQUIRE( candidates.size() == 3 );
    REQUIRE( t_type == tie_data::tie_type::K_WAY_TIE );

    candidates.clear();
    ovlp_threshold = 500;

    t_type = mod.get_tie_candidates( candidates,
                                     scores,
                                     threshold,
                                     ovlp_threshold,
                                     difference<double>()
                                   );

    REQUIRE( candidates.size() == scores.size() );
    REQUIRE( t_type == tie_data::tie_type::K_WAY_TIE );

    candidates.clear();
    threshold = 6;

    t_type = mod.get_tie_candidates( candidates,
                                     scores,
                                     threshold,
                                     ovlp_threshold,
                                     difference<double>()
                                   );

    REQUIRE( candidates.size() == scores.size() - 1 );
    REQUIRE( t_type == tie_data::tie_type::K_WAY_TIE );

    candidates.clear();

    scores.erase( scores.begin(), scores.begin() + 1 );

    ovlp_threshold = 0;

    t_type = mod.get_tie_candidates( candidates,
                                     scores,
                                     threshold,
                                     ovlp_threshold,
                                     difference<double>()
                                   );

    REQUIRE( candidates.size() == 2 );
    REQUIRE( t_type == tie_data::tie_type::TWO_WAY_TIE );

    candidates.clear();
    scores.erase( scores.begin(), scores.begin() + 1 );

    t_type = mod.get_tie_candidates( candidates,
                                     scores,
                                     threshold,
                                     ovlp_threshold,
                                     difference<double>()
                                   );

    REQUIRE( candidates.size() == 1 );
    REQUIRE( t_type == tie_data::tie_type::SINGLE_WAY_TIE );

}

TEST_CASE( "get_tie_candidates_ratio", "[module_deconv]" )
{
    auto mod = module_deconv();

    std::vector<std::pair<std::string, double>> candidates;
    std::vector<std::pair<std::string, double>> scores;
    double threshold      = 4.0;
    double ovlp_threshold = 0;

    scores.emplace_back( "100", 245.0 );
    scores.emplace_back( "101", 105.0 );
    scores.emplace_back( "102", 5.0 );
    scores.emplace_back( "103", 245.0 );
    scores.emplace_back( "104", 135.0 );
    scores.emplace_back( "105", 249.0 );

    std::sort( scores.begin(),
               scores.end(),
               util::compare_pair_non_increasing
               <std::string, double>()
            );

    auto t_type = mod.get_tie_candidates( candidates,
                                          scores,
                                          threshold,
                                          ovlp_threshold,
                                          util::ratio<double>()
                                        );

    REQUIRE( candidates.size() == 1 );
    REQUIRE( t_type == tie_data::tie_type::SINGLE_WAY_TIE );

    candidates.clear();

    ovlp_threshold = 0.5;

    t_type = mod.get_tie_candidates( candidates,
                                     scores,
                                     threshold,
                                     ovlp_threshold,
                                     util::ratio<double>()
                                   );


    REQUIRE( candidates.size() == 4 );
    REQUIRE( t_type == tie_data::tie_type::K_WAY_TIE );


    candidates.clear();

    ovlp_threshold = 0.98;

    t_type = mod.get_tie_candidates( candidates,
                                     scores,
                                     threshold,
                                     ovlp_threshold,
                                     ratio<double>()
                                   );


    REQUIRE( candidates.size() == 3 );
    REQUIRE( t_type == tie_data::tie_type::K_WAY_TIE );

    candidates.clear();

    ovlp_threshold = 0.00000001;
    threshold = 6;

    t_type = mod.get_tie_candidates( candidates,
                                     scores,
                                     threshold,
                                     ovlp_threshold,
                                     ratio<double>()
                                   );


    REQUIRE( candidates.size() == scores.size() - 1 );
    REQUIRE( t_type == tie_data::tie_type::K_WAY_TIE );

}

TEST_CASE( "test_ratio", "[struct_ratio]" )
{

    REQUIRE( ratio<double>()( 0.4, 1.0 ) == 0.4 );
    REQUIRE( ratio<double>()( 1.0, 0.4 ) == 0.4 );
    REQUIRE( ratio<double>()( 1.0, 1.0 ) == 1.0 );
    REQUIRE( ratio<double>()( 0, 1.0 ) == 0.0 );

}

TEST_CASE( "test_difference", "[struct_difference]" )
{

    REQUIRE( difference<double>()( 0.4, 1.0 ) == -0.6 );
    REQUIRE( difference<double>()( 1.0, 0.4 ) == 0.6 );
    REQUIRE( difference<double>()( 1.0, 1.0 ) == 0.0 );
    REQUIRE( difference<double>()( 0, 1.0 ) == -1.0 );

}

TEST_CASE( "calculate_overlap", "[module_deconv]" )
{

    auto mod = module_deconv();

    std::unordered_map<std::string,std::vector<std::string>> id_p_map;
    std::unordered_map<std::string,std::unordered_map<std::string,double>> counts_map;
    std::string first = "100";
    std::string second = "200";


    std::vector<std::string> f_vec { "pep1", "pep2", "pep3", "pep4" };
    std::vector<std::string> s_vec { "pep1", "pep2", "pep3", "pep4" };

    id_p_map[ first ] = f_vec;
    id_p_map[ second ] = s_vec;

    for( const auto& pep : f_vec )
        {
            counts_map[ pep ] = std::unordered_map<std::string,double>();
            counts_map[ pep ].emplace( first, 1 );
            counts_map[ pep ].emplace( second, 1 );
        }

    double threshold = 1;

    overlap_data<double> so = mod.calculate_overlap( id_p_map,
                                                     counts_map,
                                                     first, second,
                                                     evaluation_strategy::tie_eval_strategy::INTEGER_TIE_EVAL
                                                   );
    REQUIRE( so.sufficient( threshold ) );

    so = mod.calculate_overlap( id_p_map,
                                 counts_map,
                                 first, second,
                                 evaluation_strategy::tie_eval_strategy::INTEGER_TIE_EVAL
                               );

    REQUIRE( so.sufficient( 4 ) );

    so = mod.calculate_overlap( id_p_map,
                                 counts_map,
                                 first, second,
                                 evaluation_strategy::tie_eval_strategy::INTEGER_TIE_EVAL
                               );
    REQUIRE( !so.sufficient( 5 ) );


    so = mod.calculate_overlap( id_p_map,
                                 counts_map,
                                 first, second,
                                 evaluation_strategy::tie_eval_strategy::PERCENT_TIE_EVAL
                               );
    REQUIRE( so.sufficient( 0.5 ) );


}

TEST_CASE( "calculate_overlap_ratio", "[module_deconv]" )
{

    auto mod = module_deconv();

    std::unordered_map<std::string,std::vector<std::string>> id_p_map;
    std::unordered_map<std::string,std::unordered_map<std::string,double>> counts_map;
    std::string first = "100";
    std::string second = "200";


    std::vector<std::string> f_vec { "pep1", "pep3", "pep4" };
    std::vector<std::string> s_vec { "pep1", "pep2", "pep3", "pep4" };

    id_p_map[ first ] = f_vec;
    id_p_map[ second ] = s_vec;

    for( const auto& pep : f_vec )
        {
            counts_map[ pep ] = std::unordered_map<std::string,double>();
            counts_map[ pep ].emplace( first, 1 );
        }
    for( const auto& pep : s_vec )
        {
            if( counts_map.find( pep ) == counts_map.end() )
                {
                    counts_map[ pep ] = std::unordered_map<std::string,double>();
                }
            counts_map[ pep ].emplace( second, 1 );
        }



    overlap_data<double> so;

    so = mod.calculate_overlap( id_p_map,
                                counts_map,
                                first, second,
                                evaluation_strategy::tie_eval_strategy::PERCENT_TIE_EVAL
                              );

    REQUIRE( so.sufficient( 0.5 ) );

    so = mod.calculate_overlap( id_p_map,
                                counts_map,
                                first, second,
                                evaluation_strategy::tie_eval_strategy::PERCENT_TIE_EVAL
                              );
    REQUIRE( !so.sufficient( 0.8 ) );

    so = mod.calculate_overlap( id_p_map,
                                counts_map,
                                first, second,
                                evaluation_strategy::tie_eval_strategy::PERCENT_TIE_EVAL
                              );
    REQUIRE( !so.sufficient( 0.76 ) );

    so = mod.calculate_overlap( id_p_map,
                                counts_map,
                                first, second,
                                evaluation_strategy::tie_eval_strategy::PERCENT_TIE_EVAL
                              );
    REQUIRE( so.sufficient( 0.75 ) );


}

TEST_CASE( "calculate_overlap_summation", "[module_deconv]" )
{

    auto mod = module_deconv();

    std::unordered_map<std::string,std::vector<std::string>> id_p_map;
    std::unordered_map<std::string,std::unordered_map<std::string,double>> counts_map;
    std::string first = "100";
    std::string second = "200";


    std::vector<std::string> f_vec { "pep1", "pep3", "pep4" };
    std::vector<std::string> s_vec { "pep1", "pep2", "pep3", "pep4" };

    id_p_map[ first ] = f_vec;
    id_p_map[ second ] = s_vec;

    for( const auto& pep : f_vec )
        {
            counts_map[ pep ] = std::unordered_map<std::string,double>();
            counts_map[ pep ].emplace( first, 1 );
        }
    for( const auto& pep : s_vec )
        {
            if( counts_map.find( pep ) == counts_map.end() )
                {
                    counts_map[ pep ] = std::unordered_map<std::string,double>();
                }
            counts_map[ pep ].emplace( second, 1 );
        }


    overlap_data<double> so;

    so = mod.calculate_overlap( id_p_map,
                                 counts_map,
                                 first, second,
                                 evaluation_strategy::tie_eval_strategy::SUMMATION_SCORING_TIE_EVAL
                               );
    REQUIRE( so.sufficient( 0.5 ) );

    so = mod.calculate_overlap( id_p_map,
                                counts_map,
                                first, second,
                                evaluation_strategy::tie_eval_strategy::SUMMATION_SCORING_TIE_EVAL
                               );
    REQUIRE( so.sufficient( 0.75 ) );

    so = mod.calculate_overlap( id_p_map,
                                 counts_map,
                                 first, second,
                                 evaluation_strategy::tie_eval_strategy::SUMMATION_SCORING_TIE_EVAL
                               );
    REQUIRE( !so.sufficient( 0.76 ) );



}

TEST_CASE( "get_tie_type", "[module_deconv]" )
{
    auto mod = module_deconv();

    REQUIRE_THROWS( mod.get_tie_type( 0 ), std::runtime_error( "" ) );
    REQUIRE( mod.get_tie_type( 1 ) == tie_data::tie_type::SINGLE_WAY_TIE );
    REQUIRE( mod.get_tie_type( 2 ) == tie_data::tie_type::TWO_WAY_TIE );
    REQUIRE( mod.get_tie_type( 3 ) == tie_data::tie_type::K_WAY_TIE );
    REQUIRE( mod.get_tie_type( 4 ) == tie_data::tie_type::K_WAY_TIE );
}

TEST_CASE( "compare_pair", "[module_deconv]" )
{
    std::pair<int, int> p1{ 1, 4 };
    std::pair<int, int> p2{ 1, 5 };

    REQUIRE( !compare_pair_non_increasing<int,int>()( p1, p2 ) );
    REQUIRE( compare_pair_non_decreasing<int,int>()( p1, p2 ) );

    REQUIRE( compare_pair_non_increasing<int,int>()( p2, p1 ) );
    REQUIRE( !compare_pair_non_decreasing<int,int>()( p2, p1 ) );
}


TEST_CASE( "test_threshold_type", "[module_deconv]" )
{
    auto mod = module_deconv();

    REQUIRE( mod.use_ratio_score_tie_thresh( 0.5 ) );
    REQUIRE( !mod.use_ratio_score_tie_thresh( 0.0 ) );
    REQUIRE( !mod.use_ratio_score_tie_thresh( 1.0 ) );

    REQUIRE( mod.use_ratio_overlap_threshold( 0.5 ) );
    REQUIRE( !mod.use_ratio_overlap_threshold( 0.0 ) );
    REQUIRE( !mod.use_ratio_overlap_threshold( 1.0 ) );

}

TEST_CASE( "get_map_value", "[module_deconv]" )
{
    auto mod = module_deconv();

    std::unordered_map<int,int> map;
    map.emplace( 1, 1 );
    map.emplace( 2, 2 );
    map.emplace( 3, 3 );

    REQUIRE( mod.get_map_value( map, 1 ) == 1 );
    REQUIRE( mod.get_map_value( map, 2 ) == 2 );
    REQUIRE( mod.get_map_value( map, 3 ) == 3 );
    REQUIRE( mod.get_map_value( map, 4 ) == 0 );
    REQUIRE( mod.get_map_value( map, 5, 9 ) == 9 );
}

TEST_CASE( "all_distances", "[util]" )
{
    std::vector<int> data{ 1, 2, 3, 4, 5 };
    std::vector<int> distances;
    int target = 1;

    auto distance = [&]( const int a, const int b ){ return std::abs( a - b ); };

    distance_tools::all_distances( distances, data.begin(), data.end(), target, distance );

    REQUIRE( distances.size() == data.size() );

    for( std::size_t index = 0; index < distances.size(); ++index )
        {
            REQUIRE(
                    distances[ index ] == std::abs( data[ index ] - target  )
                   );
        }

}

TEST_CASE( "pairwise_distances_int", "[distance_tools]" )
{
    distance_matrix<int> dist( 5 );

    std::vector<int> data{ 1, 2, 3, 4, 5 };
    auto int_dist = []( int a, int b ){ return std::abs( a - b ); };

    distance_tools::pairwise_distances( dist,
                                        data.begin(),
                                        data.end(),
                                        int_dist
                                      );

    for( std::size_t index = 0; index < data.size(); ++index )
        {
            for( std::size_t inner_index = 0; inner_index < data.size(); ++inner_index )
                {

                    REQUIRE( dist[ index ][ inner_index ]
                             == int_dist( data[ index ], data[ inner_index ] )
                           );
                }
        }

}

TEST_CASE( "pairwise_dist_string", "[distance_tools]" )
{
    distance_matrix<int> dist( 3 );
    auto ham_dist = []( std::string a, std::string b ) -> int
        {
            std::size_t index = 0;
            int sum = 0;
            for( index = 0; index < a.length(); ++index )
                {
                    sum += a[ index ] != b[ index ];
                }
            return sum;
        };

    std::string a = "string_a";
    std::string b = "string_b";
    REQUIRE( ham_dist( a, b ) == 1 );

    std::vector<std::string> data;

    data.push_back( std::string( "string1" ) );
    data.push_back( std::string( "1234567" ) );
    data.push_back( std::string( "string3" ) );

    distance_tools::pairwise_distances( dist,
                                        data.begin(),
                                        data.end(),
                                        ham_dist
                                      );

    for( std::size_t index = 0; index < data.size(); ++index )
        {
            for( std::size_t inner_index = 0; inner_index < data.size(); ++inner_index )
                {

                    REQUIRE( dist[ index ][ inner_index ]
                             == ham_dist( data[ index ], data[ inner_index ] )
                           );
                }
        }



}

TEST_CASE( "overlap_data", "[module_deconv]" )
{
    overlap_data<double> ovlp{ 0.75, 0.5 };

    REQUIRE( ovlp.get_a_to_b() == 0.75 );
    REQUIRE( ovlp.get_b_to_a() == 0.5  );

    REQUIRE( ovlp.sufficient( 0.49 ) );
    REQUIRE( ovlp.sufficient( 0.5 ) );
    REQUIRE( !ovlp.sufficient( 0.51 ) );
    REQUIRE( ovlp.sufficient( 0 ) );


    overlap_data<int> ovlp_int{ 75, 50 };
    REQUIRE( ovlp_int.get_a_to_b() == 75 );
    REQUIRE( ovlp_int.get_b_to_a() == 50  );

    REQUIRE( ovlp_int.sufficient( 49 ) );
    REQUIRE( ovlp_int.sufficient( 50 ) );
    REQUIRE( !ovlp_int.sufficient( 51 ) );
    REQUIRE( ovlp_int.sufficient( 0 ) );

    overlap_data<int> ovlp1{ 20, 25 };
    overlap_data<int> ovlp2{ 23, 27 };

    REQUIRE(    ovlp1 < ovlp2 );
    REQUIRE( !( ovlp1 > ovlp2 ) );
    REQUIRE( !( ovlp1 == ovlp2 ) );
    REQUIRE(  ( ovlp1 != ovlp2 ) );
    REQUIRE(  ( ovlp1 <= ovlp2 ) );
    REQUIRE( !( ovlp1 >= ovlp2 ) );



}


TEST_CASE( "pair_compare", "[util]" )
{
    std::pair<int,int> pair_1{ 1, 2 };
    std::pair<int,int> pair_2{ 2, 3 };


    auto gt = []( const int& a, const int& b ){ return a > b; };
    REQUIRE( !util::pair_positional_compare( pair_1, pair_2,
                                             gt,
                                             gt
                                           )
           );

    REQUIRE( util::pair_positional_compare( pair_2, pair_1,
                                            gt,
                                            gt
                                          )
           );

    std::pair<double,std::size_t> p1{ 10.5, 11 };
    std::pair<double,std::size_t> p2{ 11.5, 10 };

    auto gt_d = []( const double& a, const double& b ){ return a > b; };
    auto gt_s = []( const std::size_t& a, const std::size_t& b ){ return a > b; };

    REQUIRE( !util::pair_positional_compare( p2, p1,
                                             gt_d,
                                             gt_s
                                           )
           );

    p1.first = 100.2;
    p1.second = 400;


    REQUIRE( !util::pair_positional_compare( p2, p1,
                                             gt_d,
                                             gt_s
                                           )
           );
    REQUIRE( util::pair_positional_compare( p1, p2,
                                             gt_d,
                                             gt_s
                                           )
           );



}

TEST_CASE( "Deconv end_to_end", "[module_deconv]" )
{
    module_deconv mod;
    options_deconv opts;

    opts.create_linkage = false;

    opts.linked_fname = std::string( "../test/test_pep_linkages.tsv" );
    opts.output_fname = std::string( "../test/test_deconv_output.tsv" );
    opts.enriched_fname = std::string( "../test/test_enriched_file.tsv" );
    opts.id_name_map_fname = std::string();

    opts.threshold = 00;
    opts.single_threaded = false;
    opts.fractional_scoring = false;
    opts.summation_scoring = true;
    opts.score_filtering = true;
    opts.species_peptides_out = "";
    opts.score_tie_threshold = 0.90;
    opts.score_overlap_threshold = 0.5;

    mod.run( &opts );
}

TEST_CASE( "empty", "[util]" )
{
    std::string iter1{ "Hello world!" };

    REQUIRE( !util::empty( iter1 ) );

    iter1 = "";

    REQUIRE( util::empty( iter1 ) );

    std::vector<int> iter2{ 1, 2, 3, 4 };
    REQUIRE( !util::empty( iter2 ) );

    iter2.clear();

    REQUIRE( util::empty( iter1 ) );

}

TEST_CASE( "global_original_scores", "[module_deconv]" )
{
    module_deconv mod;
    std::stringstream stream;

    std::unordered_map<std::string,std::string> id_name_map;
    std::unordered_map<std::string,std::string> *id_name_map_ptr = &id_name_map;


    // ordered map used so we can guarantee order
    std::map<std::string,std::pair<double,double>> score_map;

    std::vector<std::string> species{ "sp1", "sp2", "sp3", "sp4", "sp5" };

    std::string expected_include_name = "Species ID\tSpecies Name\tCount\tScore\n"
        "sp1\tsp1_name\t0\t1\n"
        "sp2\tsp2_name\t1\t2\n"
        "sp3\tsp3_name\t2\t3\n"
        "sp4\tsp4_name\t3\t4\n"
        "sp5\tsp5_name\t4\t5\n";

    std::string expected_exclude_name = "Species ID\tSpecies Name\tCount\tScore\n"
        "sp1\tsp1\t0\t1\n"
        "sp2\tsp2\t1\t2\n"
        "sp3\tsp3\t2\t3\n"
        "sp4\tsp4\t3\t4\n"
        "sp5\tsp5\t4\t5\n";


    std::size_t idx = 0;
    double idx_dbl = 1.0;
    for( auto x : species )
        {
            std::string name = x + "_name";
            id_name_map.emplace( x, name );

            score_map.emplace( x,
                               std::make_pair( idx, idx_dbl )
                             );
            ++idx;
            ++idx_dbl;
        }

    mod.write_scores( stream, id_name_map_ptr, score_map );
    REQUIRE( !stream.str().compare( expected_include_name ) );

    stream.str( std::string() );
    stream.clear();

    id_name_map_ptr = nullptr;

    mod.write_scores( stream, id_name_map_ptr, score_map );

    REQUIRE( !stream.str().compare( expected_exclude_name ) );
}


TEST_CASE( "to_dir_name", "[fs_tools]" )
{
    std::string name_no_trail = "dir1";
    std::string name_trail = "dir1/";
    std::string dest_str;

    REQUIRE( fs_tools::to_dir_name( dest_str, name_no_trail ) == name_trail );

    dest_str.clear();

    REQUIRE( fs_tools::to_dir_name( dest_str, name_trail ) == name_trail );

    dest_str.clear();
    
    REQUIRE( fs_tools::to_dir_name( name_no_trail ) == name_trail );

    dest_str.clear();

    REQUIRE( fs_tools::to_dir_name( name_trail ) == name_trail );

    std::string path_base( "dir1" );
    std::string path_base_inc( "dir1/" );

    dest_str.clear();

    fs_tools::create_fname( dest_str, path_base,
                            "fpart1", 1, 2, 3,
                            "fpart4"
                          );

    REQUIRE( dest_str == "dir1/fpart1123fpart4" );

    dest_str.clear();
    
    fs_tools::create_fname( dest_str, path_base_inc,
                            "fpart1", 1, 2, 3,
                            "fpart4"
                          );

    REQUIRE( dest_str == "dir1/fpart1123fpart4" );

}

TEST_CASE( "filter_counts (vector template)", "[module_deconv]" )
{
    std::vector<std::pair<std::size_t,double>> filter_vec; 
    module_deconv mod;

    for( std::size_t index = 0; index < 100; ++index )
        {
            filter_vec.emplace_back( std::make_pair( index, index + 1 ) );
        }


    mod.filter_counts<std::size_t,double>( filter_vec, 50  );

    REQUIRE( filter_vec.size() == 51 );
    auto comp_pair = []( const std::pair<size_t,double> first,
                         const std::pair<size_t,double> second
                         ){ return first.first < second.first; };
    REQUIRE( std::min_element( filter_vec.begin(), filter_vec.end(), comp_pair )->second == 50 );
}

TEST_CASE( "peptide", "[peptide]" )
{
    peptide pep( "pep1", "ATGC" );
    peptide pep2( "p11", "ATGC" );

    REQUIRE( pep == pep2 );

    pep2.set_sequence( "AGGG" );

    REQUIRE( pep != pep2 );
    REQUIRE( !pep2.get_sequence().compare( "AGGG" ) );

}

TEST_CASE( "scored peptide", "[peptide]" )
{
    scored_peptide<double> sc( "pep1", "ATGC", 100.0 );

    REQUIRE( sc.get_score() == 100.0 );
    REQUIRE( !sc.get_name().compare( "pep1" ) );
    REQUIRE( !sc.get_sequence().compare( "ATGC" ) );

    sc.set_score( 5 );
    REQUIRE( sc.get_score() == 5 );


}

TEST_CASE( "scored_entity", "[scored_entity.h]" )
{
    scored_entity<std::string,double> sc( std::string( "pep1" ),
                                          100.0
                                        );

    scored_entity<std::string,double> sc2( std::string( "pep1" ),
                                          155.0
                                        );

    REQUIRE( sc == sc2 );
    REQUIRE( !( sc != sc2 ) );

    REQUIRE( sc.get_score() == 100.0 );
    REQUIRE( sc2.get_score() == 155.0 );

    REQUIRE( !( sc.get_key().compare( "pep1" ) ) );
    REQUIRE( !( sc2.get_key().compare( "pep1" ) ) );

    sc.set_key( std::string( "new_key" ) );
    sc2.set_key( std::string( "new_key" ) );

    REQUIRE( sc == sc2 );
    REQUIRE( !( sc != sc2 ) );

    std::unordered_set<scored_entity<std::string,double>> se_set;

    se_set.insert( sc );

    REQUIRE( se_set.size() == 1 );

    se_set.insert( sc2 );

    REQUIRE( se_set.size() == 1 );

}

TEST_CASE( "get_kmer_frequency", "[kmer_tools]" )
{
    using scored_kmer = scored_entity<std::string,std::size_t>;
    std::unordered_set<scored_kmer> kmer_frequency_found;
    std::unordered_set<scored_kmer> kmer_frequency_expected;

    fasta_parser fp;
    std::vector<sequence> vec;
    vec = fp.parse( "../test/test.fasta" );
    std::ifstream input_stream( "../test/test_kmer_frequency.tsv",
                                std::ios_base::in
                              );


    test_utils::parse_kmer_frequency
        <std::unordered_set>( kmer_frequency_expected,
                              input_stream
                            );

    kmer_tools::get_kmer_frequencies( kmer_frequency_found,
                                      vec,
                                      9
                                    );

    REQUIRE( kmer_frequency_expected.size()
             == kmer_frequency_found.size()
           );

    for( auto current : kmer_frequency_expected )
        {
            auto found = kmer_frequency_found.find( current );
            REQUIRE( found
                     != kmer_frequency_found.end()
                   );

            REQUIRE( current.get_key()   == found->get_key() );
            REQUIRE( current.get_score() == found->get_score() );
        }


}
