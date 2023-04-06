#define CATCH_CONFIG_MAIN
#include "catch.hpp"

#include <iostream>
#include <stdexcept>
#include <fstream>
#include <algorithm>
#include <vector>
#include <unordered_map>
#include <cmath>
#include <unordered_set>
#include <string>
#include <fstream>
#include <streambuf>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>

#include "test_utils.h"
#include "overlap_data.h"
#include "distance_tools.h"
#include "sequence_indexer.h"
#include "sequence.h"
#include "fastq_sequence.h"
#include "et_search.h"
#include "maps.h"
#include "distance_matrix.h"
#include "fastq_parser.h"
#include "module_demux.h"
#include "fif_parser.h"
#include "module_deconv.h"
#include "module_info.h"
#include "module_subjoin.h"
#include "module_link.h"
#include "metadata_map.h"
#include "module_zscore.h"
#include "samplelist_parser.h"
#include "sample.h"
#include "setops.h"
#include "fastq_score.h"
#include "kmer_tools.h"
#include "fs_tools.h"
#include "peptide.h"
#include "scored_peptide.h"
#include "species_id.h"
#include "scored_entity.h"
#include "species_data.h"
#include "module_normalize.h"
#include "matrix.h"
#include "stats.h"
#include "peptide_bin.h"
#include "module_bin.h"
#include "probe_rank.h"
#include "module_enrich.h"
#include "predicate.h"
#include "file_io.h"
#include "translation_map.h"
#include "nt_aa_translator.h"

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

TEST_CASE( "fasta_parser is able to read files that exist, properly creates error when file cannot be found", "[fasta_parser]" )
{
    SECTION( "fasta_parser reads a well-formed test" )
    {
        fasta_parser fp;
        std::vector<sequence> vec;
        vec = fp.parse( "../test/input_data/test.fasta" );

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

    SECTION( "fasta_parser throws error when file is not found" )
    {
        fasta_parser fp;
        REQUIRE_THROWS( fp.parse( "does_not_exist.fasta" ) );
    }

}

TEST_CASE( "Parse Fastq", "[fastq_parser]" )
{
    std::vector<fastq_sequence> seq_vec;
    std::string fastq_fname = "../test/input_data/test_fastq.fastq";
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

TEST_CASE("Test for proper functionality of FIF Parser and Samplielist Parser", "[fif_parser], [samplelist_parser]")
{
	// define flex_idx vector
	std::vector<flex_idx> flex_idx_vec;
	// define flexible index file parser
	fif_parser fif_p;
	// parse fif file, capture result in flex_idx vector
	flex_idx_vec = fif_p.parse("../test/input_data/test_flex_index_file.tsv");

	SECTION("Flexible Index File Parser reads and creates a flex_idx vector from given file")
	{
		// verify flex_idx_vec
		auto& it = flex_idx_vec[0];
		REQUIRE(it.idx_name.compare("Index1") == 0);
		REQUIRE(it.read_name.compare("r1") == 0);
		REQUIRE(it.idx_start == 12);
		REQUIRE(it.idx_len == 10);
		REQUIRE(it.num_mismatch == 1);

		it = flex_idx_vec[1];
		REQUIRE(it.idx_name.compare("Index2") == 0);
		REQUIRE(it.read_name.compare("r2") == 0);
		REQUIRE(it.idx_start == 0);
		REQUIRE(it.idx_len == 8);
		REQUIRE(it.num_mismatch == 1);
	}

	// initialize sample list parser
	samplelist_parser slp;
	// define sample vector
	std::vector<sample> sample_vec;

	// initialize demux options
	options_demux d_opts;
    d_opts.index_fname = "../test/input_data/test_barcodes.fa";
	d_opts.samplelist_fname = "../test/input_data/test_samplelist_NS30.tsv";
	d_opts.samplename = "SampleName";	// needed for sampleliest to reference

	SECTION("Samplelist Parser can use a flex_idx vector and a properly formatted samplelist file to create a sample vector")
	{
		// parse using demux options and flex_idx vector, capture result in sample vector
		sample_vec = slp.parse(&d_opts, flex_idx_vec);

		// verify sample vector
		auto& it = sample_vec[0];
		REQUIRE(it.name.compare("Nor.27_1X_NS30_1x_ProG_A") == 0);
		REQUIRE(it.id == 0);
		REQUIRE(it.string_ids[0].compare("F_3") == 0);
		REQUIRE(it.string_ids[1].compare("R_13") == 0);

		it = sample_vec[1];
		REQUIRE(it.name.compare("Nor.28_1X_NS30_1x_ProG_A") == 0);
		REQUIRE(it.id == 1);
		REQUIRE(it.string_ids[0].compare("F_3") == 0);
		REQUIRE(it.string_ids[1].compare("R_14") == 0);

		it = sample_vec[2];
		REQUIRE(it.name.compare("BB.10_1X_NS30_1x_ProG_A") == 0);
		REQUIRE(it.id == 2);
		REQUIRE(it.string_ids[0].compare("F_11") == 0);
		REQUIRE(it.string_ids[1].compare("R_15") == 0);

		it = sample_vec[3];
		REQUIRE(it.name.compare("BB.11_1X_NS30_1x_ProG_A") == 0);
		REQUIRE(it.id == 3);
		REQUIRE(it.string_ids[0].compare("F_11") == 0);
		REQUIRE(it.string_ids[1].compare("R_16") == 0);

		it = sample_vec[4];
		REQUIRE(it.name.compare("Nor.43_1X_NS30_1x_ProG_A") == 0);
		REQUIRE(it.id == 4);
		REQUIRE(it.string_ids[0].compare("F_19") == 0);
		REQUIRE(it.string_ids[1].compare("R_13") == 0);

		it = sample_vec[5];
		REQUIRE(it.name.compare("JA.s_1X_NS30_1X_ProG_B") == 0);
		REQUIRE(it.id == 5);
		REQUIRE(it.string_ids[0].compare("F_24") == 0);
		REQUIRE(it.string_ids[1].compare("R_20") == 0);
	}
}

TEST_CASE("Full test of setops' proper functionality", "[setops]")
{
	using namespace setops;	// initialize namespace

	// testing set_intersection()
	SECTION("set_intersection() returns expected intersection")
	{
		std::set<int> int_set1 = {10, 26, 7, 9, 0};
		std::set<int> int_set2 = {10, 77, 28, 9, 100};
		std::set<int> dest_set;

		// capture result of set_intersection in dest vector
		set_intersection(dest_set, int_set1, int_set2);

		// verify dest vector is expected
		REQUIRE(dest_set.find(10) != dest_set.end());
		REQUIRE(dest_set.find(9) != dest_set.end());
	}
	/*	TODO: ask about testing this
	SECTION("set_intersection() with Get")
	{
		std::vector<int> int_vec1 = {10, 26, 7, 9, 0};
		std::vector<int> int_vec2 = {10, 77, 28, 9, 100};
		std::vector<int> dest_vec;

		set_intersection(dest_vec, int_vec1, int_vec2);

		REQUIRE(dest_vec[0] == 10);
		REQUIRE(dest_vec[1] == 9);
	}
	*/
	SECTION("set_intersection() with vectors of class K and sequential set")
	{
		std::vector<std::string> str_vec1 = {
			"pep1", "pep2", "pep3", "pep5", "pep13"
		};
		std::vector<std::string> str_vec2 = {
			"pep1", "pep10", "pep14", "pep3", "pep5"
		};
		std::vector<std::string> dest_vec;

		set_intersection(dest_vec, str_vec1, str_vec2);

		REQUIRE(dest_vec[0].compare("pep1") == 0);
		REQUIRE(dest_vec[1].compare("pep3") == 0);
		REQUIRE(dest_vec[2].compare("pep5") == 0);
	}

	// testing set_union()
	SECTION("set_union() with sequential set of type K returns expected")
	{
		std::vector<std::string> str_vec1 = {
			"pep1", "pep22", "pep5", "pep6", "pep10"
		};
		std::vector<std::string> str_vec2 = {
			"pep22", "pep6", "pep100", "pep50"
		};
		std::set<std::string> dest_set;
		std::vector<std::string> expected_set = {
			"pep1", "pep10", "pep100", "pep22", "pep5", "pep50", "pep6" 
		};

		// capture result of set_union in dest vector
		set_union(dest_set, str_vec1, str_vec2);

		// verify dest vector is expected
		std::size_t i = 0;
		for (auto dest_it : dest_set)
		{
			REQUIRE(dest_it == expected_set[i]);
			i += 1;
		}
	}
	/* TODO: ask about testing this
	SECTION("set_union() with Get")
	{
	}
	*/
	SECTION("set_union() with unordered set of type K")
	{
		std::set<int> int_set1 = {55, 6, 12, 10, 89, 100};
		std::set<int> int_set2 = {55, 16, 10, 88, 2, 0, 55};
		std::set<int> dest_set;
		std::vector<int> expected_set = {0, 2, 6, 10, 12, 16, 55, 88, 89, 100};

		set_union(dest_set, int_set1, int_set2);

		REQUIRE(dest_set.size() == expected_set.size());
		std::size_t i = 0;
		for (auto dest_it : dest_set)
		{
			REQUIRE(dest_it == expected_set[i]);
			i += 1;
		}
	}

	// testing set_difference()
	// TODO: review logic of operation
	SECTION("set_difference() returns expecetd")
	{
		std::set<int> int_set1 = {66, 12, 1, 0, 8, 16};
		std::set<int> int_set2 = {12, 16, 22, 50, 7, 66};
		std::set<int> dest_set;
		std::vector<int> expected_set = {0, 1, 8};

		// capture result of set_difference in dest vector
		set_difference(dest_set, int_set1, int_set2);

		// verify dest vector is expected
		REQUIRE(dest_set.size() == expected_set.size());
		std::size_t i = 0;
		for (auto dest_it : dest_set)
		{
			REQUIRE(dest_it == expected_set[i]);
			i += 1;
		}
	}

	// testing set_symmetric_difference()
	/*	TODO: ask if this function has been used at all
	SECTION("set_symmetric_difference() returns expected")
	{
		std::set<int> int_set1 = {13, 0, 1, 6, 11, 30, 70};
		std::set<int> int_set2 = {0, 11, 50, 60, 70, 6, 100};
		std::set<int> dest_set;

		// A - B = {13, 1, 30}
		// B - A = {50, 60, 100}
		std::vector<int> expected_set = {1, 13, 30, 50, 60, 100};

		// capture result of set_symmetric_difference()
		set_symmetric_difference<int, int>(dest_set, int_set1, int_set2);

		// verify dest vector is expected
		REQUIRE(dest_set.size() == expected_set.size());
		std::size_t i = 0;
		for (auto dest_it : dest_set)
		{
			REQUIRE(dest_it == expected_set[i]);
			i += 1;
		}
	}
	*/
}

TEST_CASE( "Add seqs to map", "[module_demux]" )
{
    fasta_parser fp;
    module_demux mod;

    std::vector<sequence> vec;
    std::size_t num_samples = 12;
    vec = fp.parse( "../test/input_data/test.fasta" );

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

TEST_CASE( "Diagnostics give a detailed count for the occurring read matches during a run", "[module_demux]" )
{
    // run demux using a real-world dataset (31250 records) with diagnostic output
    module_demux d_mod;
    options_demux d_opts;

    d_opts.input_r1_fname = std::string("../test/input_data/test_input_r1_NS30.fastq");
    d_opts.input_r2_fname = std::string("../test/input_data/test_input_r2_NS30.fastq");
    d_opts.index_fname = std::string("../test/input_data/test_barcodes.fa");
    d_opts.library_fname = std::string("../test/input_data/test_demux_library_NS30.fna");
    d_opts.samplelist_fname = std::string("../test/input_data/test_samplelist_NS30.tsv");
    d_opts.flexible_idx_fname = std::string( "" );
    d_opts.output_fname = std::string("../test/test_demux_output.tsv");
    d_opts.diagnostic_fname = std::string("../test/test_diagnostic_output.tsv");
    d_opts.read_per_loop = 80000;
    d_opts.aggregate_fname = std::string();
    d_opts.concatemer = std::string();
    d_opts.min_phred_score = 0;
    d_opts.phred_base = 33;
    d_opts.samplename = "SampleName";
    d_opts.indexes = "Index1,Index2";
    d_opts.sample_indexes = { "Index1", "Index2" };
    d_opts.set_info( &options_demux::seq_data, "41,40,2" );
    d_opts.set_info( &options_demux::index1_data, "12,10,1" );
    d_opts.set_info( &options_demux::index2_data,"0,8,1" );

    d_mod.run(&d_opts);
    // TODO: make this a section in a wider demux module test
    // check diagnostic output diagnostic file matches expected
    std::string expected = "../test/expected/test_expected_diagnostic_NS30.tsv";
    std::string actual = "../test/test_diagnostic_output.tsv";
    std::ifstream ifexpected( expected, std::ios_base::in );
    std::ifstream ifactual( actual, std::ios_base::in );
    std::string expected_line;
    std::string actual_line;

    bool lines_equal;
	while (!ifexpected.eof())
	{
		std::getline(ifexpected, expected_line);
        lines_equal = false;
        
        // TODO: find a more responsible way
	    std::ifstream ifactual(actual, std::ios_base::in);
        while (!ifactual.eof())
        {
            std::getline(ifactual, actual_line);
            if (expected_line.compare(actual_line) == 0)
            {
                lines_equal = true;
                break;
            }
        }
        ifactual.close();

        REQUIRE(lines_equal);
	}
}

// TODO: make this a section in a wider demux module test, as well
TEST_CASE("Demux output demostrates demux removes references with matching sequences", "[module_demux]")
{
	std::vector<std::string> split_line;
	std::string line;

	std::set<std::string> expected_set;
	std::set<std::string> actual_set;

	{	// collect expected file into a set
		std::ifstream ifexpected(
			"../test/expected/test_expected_demux_NS30.tsv",
			std::ios_base::in
		);

		while (!ifexpected.eof())
		{
			std::getline(ifexpected, line);
			boost::split(split_line, line, boost::is_any_of("\t"));
			expected_set.insert(split_line[0]);
		}
	}

	{	// collect actual file into a set
		std::ifstream ifactual(
			"../test/test_demux_output.tsv",
			std::ios_base::in
		);

		while (!ifactual.eof())
		{
			std::getline(ifactual, line);
			boost::split(split_line, line, boost::is_any_of("\t"));
			actual_set.insert(split_line[0]);
		}
	}

	// compare expected and actual sets
	REQUIRE(actual_set.size() == expected_set.size());

	auto expected_ref = expected_set.find("NS30_000000-1");
	REQUIRE(actual_set.find(*expected_ref) != actual_set.end());

	expected_ref = expected_set.find("NS30_000000-2");
	REQUIRE(actual_set.find(*expected_ref) != actual_set.end());
	
	expected_ref = expected_set.find("NS30_000001-2");
	REQUIRE(actual_set.find(*expected_ref) != actual_set.end());
	
	expected_ref = expected_set.find("NS30_000001-1");
	REQUIRE(actual_set.find(*expected_ref) != actual_set.end());
	
	expected_ref = expected_set.find("NS30_000000-3");
	REQUIRE(actual_set.find(*expected_ref) != actual_set.end());
}

// TODO: figure out how to capture stdout so info module logs do not appear
// while running pepsirf_test
/*
TEST_CASE("Full test of info module", "[module_info]")
{
	// initialize info module
    module_info i_mod;
    options_info i_opts;

    i_opts.in_fname = "../test/input_data/test_info_score_matrix.tsv";
    i_opts.out_samples_fname = "../test/test_info_sample_names.tsv";
    i_opts.out_pep_names_fname = "../test/test_info_pep_names.tsv";
    i_opts.out_col_sums_fname = "../test/test_info_col_sums.tsv";
    i_opts.in_replicates_fname = "../test/input_data/test_info_rep_names.tsv";
    i_opts.out_avgs_fname = "../test/test_info_avg_matrix.tsv";

    // run info module
    i_mod.run(&i_opts);
	
    std::string expected_line = "";
    std::string actual_line = "";
	SECTION("info module creates a properly formatted sample name file")
	{
        std::ifstream ifexpected(
            "../test/expected/test_expected_info_sample_names.tsv",
            std::ios_base::in
        );
        std::ifstream ifactual(
            "../test/test_info_sample_names.tsv",
            std::ios_base::in
        );

        while (!ifexpected.eof() && !ifactual.eof())
        {
            std::getline(ifexpected, expected_line);
            std::getline(ifactual, actual_line);
            REQUIRE(expected_line.compare(actual_line) == 0);
        }
	}
	SECTION("info module creates a properly formatted peptide name file")
	{
        std::ifstream ifexpected(
            "../test/expected/test_expected_info_pep_names.tsv",
            std::ios_base::in
        );
        std::ifstream ifactual(
            "../test/test_info_pep_names.tsv",
            std::ios_base::in
        );

        while (!ifexpected.eof() && !ifactual.eof())
        {
            std::getline(ifexpected, expected_line);
            std::getline(ifactual, actual_line);
            REQUIRE(expected_line.compare(actual_line) == 0);
        }
	}
    SECTION("info module properly calculates and formats column sums file")
    {
        std::ifstream ifexpected(
            "../test/expected/test_expected_info_col_sums.tsv",
            std::ios_base::in
        );
        std::ifstream ifactual(
            "../test/test_info_col_sums.tsv",
            std::ios_base::in
        );

        while (!ifexpected.eof() && !ifactual.eof())
        {
            std::getline(ifexpected, expected_line);
            std::getline(ifactual, actual_line);
            REQUIRE(expected_line.compare(actual_line) == 0);
        }
    }
	SECTION("info module creates and calculates average matrix file")
	{
        std::ifstream ifexpected(
            "../test/expected/test_expected_info_avg_matrix.tsv",
            std::ios_base::in
        );
        std::ifstream ifactual(
            "../test/test_info_avg_matrix.tsv",
            std::ios_base::in
        );

        while (!ifexpected.eof() && !ifactual.eof())
        {
            std::getline(ifexpected, expected_line);
            std::getline(ifactual, actual_line);
            REQUIRE(expected_line.compare(actual_line) == 0);
        }
	}
}
*/

TEST_CASE( "samplelist_parser is able to read files that exist, properly creates errors when file cannot be found/read", "[samplelist_parser]" )
{
    // SECTION( "samplelist_parser is able to read a well-formed test")
    // {
    //     samplelist_parser slp;
    //     std::string filename = "../test/test_samplelist.tsv";

    //     auto vec = slp.parse( filename );

    //     REQUIRE( vec.size() == 96 );


    //     unsigned int index = 0;
    //     std::unordered_map<sample, int> s_map( vec.size() );
    //     for( index = 0 ; index < vec.size(); ++index )
    //         {
    //             s_map[ vec[ index ] ] = vec[ index ].id;
    //         }

    //     REQUIRE( s_map.size() <= vec.size() );
    // }


    // SECTION( "samplelist_parser throws an error when a file is not found or incorrect format" )
    // {
    //     samplelist_parser sl;
    //     REQUIRE_THROWS( sl.parse( "../test/test.fasta" ) );
    //     REQUIRE_THROWS( sl.parse( "does_not_exist.tsv" ) );
    // }
}

TEST_CASE( "Test String Indexing", "[string_indexer]" )
{
    std::string origin = "";
    sequence_indexer si;
    fastq_parser fp;
    auto seqs = fp.parse( "../test/input_data/test_fastq.fastq" );
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

TEST_CASE( "Reference-independent Demultiplexing" )
{
    std::vector<sequence>
        values
    {
        sequence{ "", "AAAA" },
        sequence{ "", "AAAT" },
        sequence{ "", "AAAG" },
        sequence{ "", "AAAC" }
    };

    sequence_indexer si;

    si.index( values );

    SECTION( "Value-based value-items in the map" )
        {
            using mapped_t = std::vector<std::size_t>;
            using map_t = std::unordered_map<sequence,
                                             mapped_t
                                            >;

            map_t map;

            std::for_each( values.begin(),
                           values.end() - 1,
                           [&]( const sequence& seq )
                           {
                               mapped_t val{ 5, 0 };
                               map.emplace( seq, val );
                           }
                         );
            et_seq_search<map_t, false> search{ si, map, 2 };
            auto iter = search.find( sequence( "", "AAAA" ), 4, 0, 4 );

            REQUIRE( iter != map.end() );
            iter = search.find( sequence( "", "BAAA" ), 4, 0, 4 );
            REQUIRE( iter != map.end() );
            REQUIRE( map.size() == 4 );

            auto it = map[ sequence( "", "BAAA" ) ];
            REQUIRE( it.size() == 2 );
            REQUIRE( it[ 0 ] == 0 );
            REQUIRE( it[ 1 ] == 0 );



        }

    SECTION( "Pointer-based value-items in the map" )
        {
            using mapped_t = std::vector<std::size_t>*;
            using map_t = std::unordered_map<sequence,
                                             mapped_t
                                            >;

            map_t map;

            std::for_each( values.begin(),
                           values.end() - 1,
                           [&]( const sequence& seq )
                           {
                               mapped_t val = new std::vector<std::size_t>();
                               val->resize( 2 );
                               val->at( 0 ) = 5;
                               val->at( 1 ) = 0;
                               map.emplace( seq, val );
                           }
                         );
            et_seq_search<map_t, false> search{ si, map, 2 };
            auto iter = search.find( sequence( "", "AAAA" ), 4, 0, 4 );

            REQUIRE( iter != map.end() );
            iter = search.find( sequence( "", "BAAA" ), 4, 0, 4 );
            REQUIRE( iter != map.end() );
            REQUIRE( map.size() == 4 );

            auto it = map[ sequence( "", "BAAA" ) ];
            REQUIRE( it->size() == 2 );
            REQUIRE( it->at( 0 ) == 0 );
            REQUIRE( it->at( 1 ) == 0 );

        }

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

    vec_a = fp.parse( "../test/input_data/test_fastq.fastq" );

    std::for_each( vec_a.begin(), vec_a.end(),
                   [&]( const sequence& seq )
                   {
                       vec.emplace_back( seq.name, seq.seq );
                   }
                 );


    lib_idx.index( vec );

    parallel_map<sequence, std::vector<std::size_t>*> my_map;
    using map_t = decltype( my_map );

    mod.add_seqs_to_map( my_map, vec, num_samples );

    et_seq_search<map_t>
        searcher{ lib_idx, my_map, num_samples };

    parallel_map<sequence, std::vector<std::size_t>*>::iterator seq_match;

    std::size_t index = 0;

    // perfect matches
    for( index = 0; index < my_map.size(); ++index )
        {
            sequence a = sequence( "", vec[ index ].seq );

            seq_match = searcher.find( a,
                                       num_mismatches,
                                       seq_start,
                                       seq_length
                                     );
            REQUIRE( seq_match != my_map.end() );
            REQUIRE( seq_match->second->size() == 1 );

            seq_match->second->at( 0 )++;

            REQUIRE( seq_match->second->at( 0 ) == 1 );

        }
    sequence seq =
    sequence( "", "AGCCAGCTTGCGGCAAAACTGCGTAACCGTCTTCTCGTTCTCTAAAAACCATTTTTCGTCCCCTTCGGGGCGGTGGTCTATAGTGTTATTAATATCAAGTTGGGGGAGCACATTGTAGCATTGTGCCAATTCATCCATTAACTTCTCAGT" );

    seq_match = searcher.find( seq,
                               num_mismatches,
                               seq_start,
                               seq_length
                             );

    REQUIRE( seq_match == my_map.end() );

    seq_match = searcher.find( seq,
                               1,
                               seq_start,
                               seq_length
                             );

    REQUIRE( seq_match != my_map.end() );
    seq =
    sequence( "", "AGCCAGCTTGCGGCAAAACTGCGTAACCGTCTTCTCGTTCTCTAAAAACCATTTTTCGTCCCCTTCGGGGCGGTGGTCTATAGTGTTATTAATATCAAGTTGGGGGAGCACATTGTAGCATTGTGCCAATTCATCCATTAACTTCGCAAT" );

    seq_match = searcher.find( seq,
                               3,
                               seq_start,
                               seq_length
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

    opts.linked_fname = std::string( "../test/input_data/test_pep_linkages.tsv" );
    opts.output_fname = std::string( "../test/test_deconv_output.tsv" );
    opts.enriched_fname = std::string( "../test/input_data/test_enriched_file.tsv" );
    opts.id_name_map_fname = std::string();

    opts.threshold = 00;
    opts.single_threaded = false;
    opts.scoring_strategy = "summation";
    opts.score_filtering = true;
    opts.species_peptides_out = "";
    opts.score_tie_threshold = 0.90;
    opts.score_overlap_threshold = 0.5;

    mod.run( &opts );
}

TEST_CASE("Deconv module sorts ties alphabetically", "[module_deconv]")
{
	// initialize deconv components
	
	// run deconv module
	
	// check created files against expected files -- that a lot of data
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
    scored_entity<std::string, double> sc( std::string( "pep1" ),
                                          100.0
                                        );

    scored_entity<std::string, double> sc2( std::string( "pep1" ),
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

    std::unordered_set<scored_entity<std::string, double>> se_set;

    se_set.insert( sc );

    REQUIRE( se_set.size() == 1 );

    se_set.insert( sc2 );

    REQUIRE( se_set.size() == 1 );

	sc.set_score(226.0);
	sc2.set_score(77.0);

	REQUIRE(sc.get_score() == 226.0);
	REQUIRE(sc2.get_score() == 77.0);
}

TEST_CASE( "get_kmer_frequency", "[kmer_tools]" )
{
    using scored_kmer = scored_entity<std::string,std::size_t>;
    std::unordered_set<scored_kmer> kmer_frequency_found;
    std::unordered_set<scored_kmer> kmer_frequency_expected;

    fasta_parser fp;
    std::vector<sequence> vec;
    vec = fp.parse( "../test/input_data/test.fasta" );
    std::ifstream input_stream( "../test/input_data/test_kmer_frequency.tsv",
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

TEST_CASE( "species_data", "[module_deconv]" )
{
    species_id<std::string> spec( "species 1" );
    scored_peptide<double> peptide( "ATGC", 105.0 );
    double score = 4;
    double count = 11;

    species_data dat( spec, score, count, peptide );

    auto pep = dat.get_highest_scoring_peptide();

    REQUIRE( pep.get_score() == 105.0 );
    REQUIRE( !pep.get_sequence().compare( "ATGC" ) );
    dat.set_score( 104.0 );

    REQUIRE( dat.get_score() == 104.0 );

}

TEST_CASE( "score_peptide_for_species", "[module_deconv]" )
{
    auto mod = module_deconv();

    peptide pep( "ATGC" );

    std::unordered_map<
        std::string,std::vector<std::pair<std::string,double>>>
        spec_count_map;

    std::string species_id = "spec 1";

    evaluation_strategy::score_strategy strat;

    strat = evaluation_strategy::score_strategy::INTEGER_SCORING;

    std::vector<std::pair<std::string,double>> value;
    value.emplace_back( species_id, 150.0 );
    value.emplace_back( "spec 2", 140.0 );

    spec_count_map[ pep.get_sequence() ] = value;

    auto eval = [&]() -> double
    {
        return mod.score_peptide_for_species( pep, spec_count_map,
                                              species_id,
                                              strat
                                            );
    };
    double score = eval();

    REQUIRE( score == 1 );

    strat = evaluation_strategy::score_strategy::SUMMATION_SCORING;

    score = eval();

    REQUIRE( score == 150.0 );

    strat = evaluation_strategy::score_strategy::FRACTIONAL_SCORING;

    score = eval();

    REQUIRE( score == 0.5 );

}

TEST_CASE( "score_species_peptides/get_highest_score_per_species", "[module_deconv]" )
{
    auto mod = module_deconv();

    std::unordered_map<std::string,
                       std::vector<scored_peptide<double>>
                       >
        counts;


    std::unordered_map<std::string,std::vector<std::string>>
        id_count_map;

    std::unordered_map<std::string,std::vector<std::pair<std::string,double>>>
        spec_count_map;

     evaluation_strategy::score_strategy strat;

     std::vector<std::string> a_pep;
     std::vector<std::string> b_pep;

     a_pep.emplace_back( "pep1" );
     a_pep.emplace_back( "pep2" );

     b_pep.emplace_back( "pep1" );
     b_pep.emplace_back( "pep2" );

     id_count_map.emplace(
                          std::make_pair(
                                         "species 1",
                                         a_pep
                                         )
                          );
     id_count_map.emplace(
                          std::make_pair(
                                         "species 2",
                                         b_pep
                                         )
                          );

     std::vector<std::pair<std::string,double>> a_score;
     std::vector<std::pair<std::string,double>> b_score;

     a_score.emplace_back( std::make_pair( "species 1", 2.4 ) );
     a_score.emplace_back( std::make_pair( "species 2", 3.4 ) );

     b_score.emplace_back( std::make_pair( "species 1", 8.4 ) );
     b_score.emplace_back( std::make_pair( "species 2", 9.4 ) );

     spec_count_map[ "pep1" ] = a_score;
     spec_count_map[ "pep2" ] = b_score;

     std::unordered_map<std::string,scored_peptide<double>> highest_scores;

     auto eval_sc_sp_pep = [&]()
         {
             mod.score_species_peptides( counts,
                                         id_count_map,
                                         spec_count_map,
                                         strat
                                         );
         };

     auto clear = [&]() { counts[ "species 1" ].clear();
                          counts[ "species 2" ].clear();
                          highest_scores.clear();
     };

     auto eval_max_score = [&] ()
         {
             mod.get_highest_score_per_species( highest_scores, counts );
         };

     strat = evaluation_strategy::score_strategy::SUMMATION_SCORING;

     eval_sc_sp_pep();

     REQUIRE( counts[ "species 1" ][ 0 ].get_score() == 2.4 );
     REQUIRE( counts[ "species 1" ][ 1 ].get_score() == 8.4 );
     REQUIRE( counts[ "species 2" ][ 0 ].get_score() == 3.4 );
     REQUIRE( counts[ "species 2" ][ 1 ].get_score() == 9.4 );

     eval_max_score();

     REQUIRE( highest_scores[ "species 1" ].get_score() == 8.4 );
     REQUIRE( highest_scores[ "species 2" ].get_score() == 9.4 );

     clear();
     strat = evaluation_strategy::score_strategy::INTEGER_SCORING;

     eval_sc_sp_pep();

     REQUIRE( counts[ "species 1" ][ 0 ].get_score() == 1 );
     REQUIRE( counts[ "species 1" ][ 1 ].get_score() == 1 );
     REQUIRE( counts[ "species 2" ][ 0 ].get_score() == 1 );
     REQUIRE( counts[ "species 2" ][ 1 ].get_score() == 1 );

     eval_max_score();

     REQUIRE( highest_scores[ "species 1" ].get_score() == 1 );
     REQUIRE( highest_scores[ "species 2" ].get_score() == 1 );

     clear();
     strat = evaluation_strategy::score_strategy::FRACTIONAL_SCORING;

     eval_sc_sp_pep();

     REQUIRE( counts[ "species 1" ][ 0 ].get_score() == 0.5 );
     REQUIRE( counts[ "species 1" ][ 1 ].get_score() == 0.5 );
     REQUIRE( counts[ "species 2" ][ 0 ].get_score() == 0.5 );
     REQUIRE( counts[ "species 2" ][ 1 ].get_score() == 0.5 );

     eval_max_score();
     REQUIRE( highest_scores[ "species 1" ].get_score() == 0.5 );
     REQUIRE( highest_scores[ "species 2" ].get_score() == 0.5 );
}

TEST_CASE("Full test of normalize module", "[module_normalize]")
{
	module_normalize norm_mod;

	SECTION("Test get_sum() works as expected")
	{
		// define vec of doubles
		// initialize matrix of doubles

		// get sums from matrix

		// check destination vector for expected values
	}
	SECTION("Test get_neg_average() works as expected")
	{
		// initialize peptide data sample major
		// initialize negative filter
		// define destination map for averages

		// get negative averages

		// check destination map contains expected averages
	}
	SECTION("Test constant_factor_normalization() works as expected")
	{
		SECTION("Normalize column sums by 2")
		{
			// initialize column sums vector

			// normalize column vec by 2

			// check column vec sums have been normalized as expected
		}
		SECTION("Normalize column sums by 4")
		{
			// initialize column sums vector

			// normalize column vec by 4

			// check column vec sums have been normalized as expected
		}
		SECTION("Normalize column sums by 100")
		{
			// initialize column sums vector

			// normalize column vec by 100

			// check column vec sums have been normalized as expected
		}
	}
	SECTION("Test filter_neg_control_start() works as expected")
	{
		// TODO: filter_neg_control_start() needs a description block
	}
	SECTION("Test normalize_counts() works as expected")
	{	// TODO: should we have more of this test?
		// initialize matrix of origianl scores
		// initialize normalize factors vector (double)

		// normalize original scores

		// check original scores have been normalized as expected
	}
	SECTION("Test compute_size_factors() works as expected")
	{
		// define destination vector for size factors
		// initialize matrix of counts

		// compute size factors of data

		// check destination vec for expected size factors
	}
	SECTION("Test compute_diff() works as expected")
	{
		// define destination vector for norm score differences
		// initialize map for means of negative controls

		// compute differences

		// check destination vec for expected differences
	}
	SECTION("Test compute_diff_ratio() works as expected")
	{
		// define destination peptide data sample major
		// initialize map for means of provided controls

		// compute difference ratio

		// check destination for expected difference ratios
	}
	SECTION("Test compute_ratio() works as expected")
	{
		// define destination peptide data sample major
		// initialize map for means of provided controls

		// compute ratio

		// check destination for expected ratios
	}
}

TEST_CASE( "geometric means", "[stats]" )
{
    // TODO: RE-implement this test
    module_normalize mod;

    matrix<double> values{ 3, 4 };

    for( int outer = 0; outer < 3; ++outer )
        {
            for( int inner = 0; inner < 4; ++inner )
                {
                    values( outer, inner ) = outer * inner;
                }
        }

    double mean = stats::geom_mean( values.row_begin( 1 ),
                                    values.row_end( 1 )
                                  );
    double epsilon = 0.0005;

    // check that we're sufficiently close, I calculated the
    // value separately, taking care to remove any zeros
    REQUIRE( std::abs( mean - 1.81712 ) < epsilon );


    mean = stats::geom_mean( values.row_begin( 2 ),
                             values.row_end( 2 )
                           );

    REQUIRE( std::abs( mean - 3.63424 ) < epsilon );

    mean = stats::geom_mean( values.row_begin( 0 ),
                             values.row_end( 0 )
                           );

    REQUIRE( mean == 0 );
}

TEST_CASE( "matrix creation, setting individual members of matrix", "[matrix]" )
{
    matrix<int> a{ 12, 13 };

    a.at( 0, 0 ) = 0;

    REQUIRE( a.at( 0, 0 ) == 0 );
    REQUIRE( a( 0, 0 ) == 0 );
    REQUIRE( a.get_shape() == std::pair<std::uint32_t,std::uint32_t>( 12, 13 ) );
    REQUIRE_THROWS( a( 1300, 1400 ) );

    a( 0, 3 ) = 4;
    REQUIRE( a.at( 0, 3 ) == 4 );
    REQUIRE( a( 0, 3 ) == 4 );

    a.set_all( 100 );

    for( std::uint32_t x = 0; x < 12; ++x )
        {
            for( std::uint32_t y = 0; y < 13; ++y )
                {
                    REQUIRE( a( x, y ) == 100 );
                }
        }

    a.at( 1, 2 ) = 4;
    REQUIRE( a( 1, 2 ) == 4 );
}

TEST_CASE( "labeled_matrix", "[labeled_matrix]" )
{
    std::vector<std::string> row_labels{ "row_1", "row_2", "row_3", "row_4", "row_5" };
    std::vector<std::string> col_labels{ "col_1", "col_2", "col_3", "col_4", "col_5" };

    labeled_matrix<double, std::string> lab_mat( 5, 5, row_labels, col_labels );

    lab_mat.set_all( 15 );

    REQUIRE( lab_mat.get_shape() == std::pair<std::uint32_t,std::uint32_t>( 5, 5 ) );
    REQUIRE( lab_mat( "row_1", "col_1" ) == 15 );

    lab_mat( "row_1", "col_3" ) = 4;
    REQUIRE( lab_mat( 0, 2 ) == 4 );

    REQUIRE_THROWS( lab_mat.at( "nonsense", "labels" ) );
    REQUIRE_THROWS( lab_mat.at( 150, 2343 ) );

    std::vector<std::string> rows{ "row_1" };
    auto new_matr = lab_mat.filter_rows( rows );

    REQUIRE( new_matr.get_shape() == std::pair<std::uint32_t,std::uint32_t>( 1, 5 ) );

    for( std::uint32_t col_idx = 0; col_idx < col_labels.size(); ++col_idx )
        {
            REQUIRE( new_matr( 0, col_idx ) == lab_mat( 0, col_idx ) );
        }

    std::vector<std::string> rows2{ "row_1", "row_3" };
    auto new_matr2 = lab_mat.filter_rows( rows2 );

    REQUIRE( new_matr2.get_shape() == std::pair<std::uint32_t,std::uint32_t>( 2, 5 ) );

    for( std::uint32_t col_idx = 0; col_idx < col_labels.size(); ++col_idx )
        {
            REQUIRE( new_matr2( 0, col_idx ) == lab_mat( 0, col_idx ) );
            REQUIRE( new_matr2( 1, col_idx ) == lab_mat( 2, col_idx ) );
        }

    for( uint in_idx = 0; in_idx < 5; ++in_idx )
        {
            for( uint out_idx = 0; out_idx < 5; ++out_idx )
                {
                    lab_mat( in_idx, out_idx ) = in_idx * out_idx;
                }
        }
    std::vector<std::string> rows3{ "row_1", "row_3", "row_5", "row_4" };
    auto new_matr3 = lab_mat.filter_rows( rows3 );
    new_matr.compare_row( lab_mat, 1 );

    for( std::uint32_t col_idx = 0; col_idx < col_labels.size(); ++col_idx )
        {
            REQUIRE( new_matr3( 0, col_idx ) == lab_mat( 0, col_idx ) );
            REQUIRE( new_matr3( 1, col_idx ) == lab_mat( 2, col_idx ) );
            REQUIRE( new_matr3( 2, col_idx ) == lab_mat( 4, col_idx ) );
            REQUIRE( new_matr3( 3, col_idx ) == lab_mat( 3, col_idx ) );
        }
}

TEST_CASE( "labeled_matrix full outer join", "[matrix]" )
{
    std::uint32_t a_N = 3, a_M = 3;
    std::uint32_t b_N = 2, b_M = 6;

    std::vector<std::string> a_rl{ "arow1", "arow2", "arow3" };
    std::vector<std::string> b_rl{ "brow1", "brow2" };

    std::vector<std::string> a_cl{ "acol1", "acol2", "acol3" };
    std::vector<std::string> b_cl{ "bcol1", "bcol2", "bcol3",
                                   "bcol4", "bcol5", "bcol6"
                                 };


    labeled_matrix<double, std::string> a( a_N, a_M,
                                           a_rl, a_cl
                                         );
    labeled_matrix<double, std::string> b( b_N, b_M,
                                           b_rl, b_cl
                                         );
    a( "arow1", "acol1" ) = 5.5;
    b( "brow2", "bcol5" ) = 5.5;

    labeled_matrix<double,std::string> a_join_b
        = a.full_outer_join( b );

    REQUIRE( a_join_b.nrows() == a_N + b_N );
    REQUIRE( a_join_b.ncols() == a_M + b_M );
    REQUIRE( a_join_b( "arow1", "acol1" ) == 5.5 );
    REQUIRE( a_join_b( "brow2", "bcol5" ) == 5.5 );
    REQUIRE( a_join_b( "brow2", "acol1" ) == 0 );

    std::vector<std::string> a1_rl;
    std::vector<std::string> a1_cl;
    std::vector<std::string> a2_rl;
    std::vector<std::string> a2_cl;
    int x_dim = 100;
    int y_dim = 150;

    labeled_matrix<double,std::string> a1_matr( x_dim, y_dim );
    labeled_matrix<double,std::string> a2_matr( x_dim, y_dim );

    for( int x = 0; x < x_dim - 1; ++x )
        {

            a1_rl.push_back( "a_" );
            a1_rl.back().append( std::to_string( x ) );

            a2_rl.push_back( "c_" );
            a2_rl.back().append( std::to_string( x ) );
        }

    a1_rl.push_back( "shared_label_row" );
    a2_rl.push_back( "shared_label_row" );

    for( int y = 0; y < y_dim; ++y )
        {

            a1_cl.push_back( "b_" );
            a1_cl.back().append( std::to_string( y ) );


            a2_cl.push_back( "d_" );
            a2_cl.back().append( std::to_string( y ) );


        }

    for( int x = 0; x < x_dim; ++x )
        {
            for( int y = 0; y < y_dim; ++ y )
                {
                    a1_matr( x, y ) = x * y;
                    a2_matr( y % x_dim, x ) = x * y;
                }
        }

    a1_matr.set_col_labels( a1_cl );
    a1_matr.set_row_labels( a1_rl );
    a2_matr.set_col_labels( a2_cl );
    a2_matr.set_row_labels( a2_rl );

    labeled_matrix<double,std::string> a1_join_b2 = a1_matr.full_outer_join( a2_matr );

    for( int x = 0; x < x_dim; ++x )
        {
            for( int y = 0; y < y_dim; ++y )
                {
                    REQUIRE( a1_join_b2( a1_rl[ x ], a1_cl[ y ] ) == a1_matr( a1_rl[ x ], a1_cl[ y ] ) );
                    REQUIRE( a1_join_b2( a2_rl[ x ], a2_cl[ y ] ) == a2_matr( a2_rl[ x ], a2_cl[ y ] ) );
                }
        }


}

TEST_CASE("Testing in_range(), swap(), and filter_cols() funcitons of a labeled matrix", "labeled_matrix")
{
	// initialize row and column labels vec
	// define a 3x3 labeled matrix with the row and column label vecs
	
	// fill matrix with numbers
	
	// check element (1, 1) is in matrix
	// check element (4, 4) is in matrix
	// check element (0, 2) is in matrix
	// check element (5, 0) is in matrix
	
	// initialize new row and column labels vec
	// define new 3x3 labeled matrix
	
	// fill new matrix with zeros
	
	// swap matrices
	
	// define empty destination vector of strings
	// check original matrix row labels are the new row labels
	// clear dest vec

	// check original matrix column labels are the new column labels
	// clear dest vec

	// check new matrix row labels are the old row labels
	// clear dest vec

	// check new matrix column labels are the old column lables

	// check old matrix values are all zero
	// check new matrix values are original values
}

TEST_CASE( "median", "[stats]" )
{

    SECTION( "Median can be found in a vector of odd length" )
        {
            std::vector<int> data{ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
            auto median = stats::median( data.begin(), data.end() );
            REQUIRE( median == 5 );

            data = std::vector<int>{ 4, 5, 0, 1, 3, 2, 8, 10, 9, 7, 6 };
            median = stats::median( data.begin(), data.end() );
            REQUIRE( median == 5 );

            data = std::vector<int>{ 6, 8, 10, 2, 3, 5, 4, 7, 0, 1, 9 };
            median = stats::median( data.begin(), data.end() );
            REQUIRE( median == 5 );

            // data generated with python
            data = std::vector<int>{ 67,69,94,77,49,89,86,61,23,87,75,80,78,81,1,93,27,22,88,0,11 };
            median = stats::median( data.begin(), data.end() );
            REQUIRE( median == 75 );

        }

    SECTION( "Median can be found in a vector of even length" )
        {
            std::vector<int> data{ 0, 1, 2, 3, 4, 6, 7, 8, 9, 10 };
            auto median = stats::median( data.begin(), data.end() );
            REQUIRE( median == 5 );

            data = std::vector<int>{ 4, 5, 0, 1, 3, 2, 10, 9, 7, 6 };
            median = stats::median( data.begin(), data.end() );
            REQUIRE( median == 4.5 );

            data = std::vector<int>{ 6, 8, 10, 2, 3, 5, 4, 7, 1, 9 };
            median = stats::median( data.begin(), data.end() );
            REQUIRE( median == 5.5 );

            // data generated with python
            data = std::vector<int>{ 67,69,94,77,49,89,86,61,23,87,75,80,78,81,1,93,27,22,88,0 };
            median = stats::median( data.begin(), data.end() );
            REQUIRE( median == 76 );
        }

}

TEST_CASE( "Arithmetic mean", "[stats]" )
{

    SECTION( "Arithmetic mean of vector" )
        {
            std::vector<double> data{ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
            REQUIRE( stats::arith_mean( data.begin(), data.end() ) == 5.5 );
        }

}

TEST_CASE( "matrix transposition", "[matrix]" )
{
    SECTION( "Transposition of a square matrix", "[matrix]" )
        {
            labeled_matrix<int,std::string> matr{ 5, 5 };
            matr.set_col_labels( std::vector<std::string>{ "col_1", "col_2", "col_3", "col_4", "col_5" } );
            matr.set_row_labels( std::vector<std::string>{ "row_1", "row_2", "row_3", "row_4", "row_5" } );
            matr.set_all( 5 );

            matr( 1, 4 ) = 12;
            auto matr_t = matr.transpose();
            auto matr_t_t = matr_t.transpose();

            REQUIRE( matr_t( 4, 1 ) == 12 );

            REQUIRE( ( matr_t_t == matr ) );
        }

    SECTION( "Transposition of a rectangular matrix", "[matrix]" )
        {
            // a non-symmetric function of two variables, x and y,
            // used because for symmetric F F(x,y) = F(y,x), making it not
            // useful for testing our transposition operations
            auto non_symmetric = []( int x, int y ) -> int
                {
                    constexpr int A = 4, B = 8, R = 28;
                    return ( A * x * x ) + ( B * y * y ) - ( R * R );
                };

            labeled_matrix<int,std::string> matr{ 8, 2 };
            matr.set_row_labels( std::vector<std::string>{ "row_1", "row_2", "row_3", "row_4",
                                                           "row_5", "row_6", "row_7", "row_8"
                                                         }
                               );

            matr.set_col_labels( std::vector<std::string>{ "col_1", "col_2" } );

            for( int i = 0; i < 8; ++i )
                {
                    for( int j = 0; j < 2; ++j )
                        {
                            matr( i, j ) = non_symmetric( i, j );
                        }
                }

            auto matr_t = matr.transpose();
            auto matr_t_t = matr_t.transpose();

            for( int i = 0; i < 8; ++i )
                {
                    for( int j = 0; j < 2; ++j )
                        {
                            REQUIRE( matr_t( j, i ) == matr( i, j ) );
                        }
                }

            REQUIRE( ( matr_t_t == matr ) );
        }

}

TEST_CASE( "Standard deviation", "[stats]" )
{
    std::vector<double> data{ 10, 12, 23, 23, 16, 23, 21, 16 };

    double stdev = stats::stdev( data.begin(), data.end() );
    REQUIRE( ( stdev && stats::is_close( stdev, 5.237 ) ) );
}

TEST_CASE( "z-scores", "[stats]" )
{
    std::vector<double> data{ 10, 12, 23, 23, 16, 23, 21, 16 };

    auto verify_vec = [&]( const std::vector<double>& stdev )
        {
            REQUIRE( stats::is_close( stdev[ 0 ], -1.528 ) );
            REQUIRE( stats::is_close( stdev[ 1 ], -1.1456 ) );
            REQUIRE( stats::is_close( stdev[ 2 ], 0.9547  ) );
            REQUIRE( stats::is_close( stdev[ 3 ], 0.9547  ) );
            REQUIRE( stats::is_close( stdev[ 4 ], -0.38189 ) );
            REQUIRE( stats::is_close( stdev[ 5 ], 0.9547  ) );
            REQUIRE( stats::is_close( stdev[ 6 ], 0.5728  ) );
            REQUIRE( stats::is_close( stdev[ 7 ], -0.38189 ) );
        };

    SECTION( "Input and output are not the same vector" )
        {
            std::vector<double> deviations;
            deviations.resize( 8 );

            stats::zscore( data.begin(), data.end(), deviations.begin() );
            verify_vec( deviations );
        }

    SECTION( "Input and output are the same vector" )
        {
            stats::zscore( data.begin(), data.end(), data.begin() );
            verify_vec( data );
        }
}

TEST_CASE( "Parsing/Writing Bins from stream", "[peptide_bin]" )
{
    std::stringstream bins_in;
    bins_in << "pep_1\tpep_2\tpep_3\npep_4\tpep_5\tpep_6\n";
    bin_collection bin_c = peptide_bin_io::parse_bins( bins_in );

    std::stringstream bins_out;

    peptide_bin_io::write_bins( bins_out, bin_c );

	// test parse_bins() -> tests add_bin()
    REQUIRE( ( bin_c == peptide_bin_io::parse_bins( bins_out ) ) );

	// test contains() method
    auto first_bin = bin_c.begin();
    REQUIRE( first_bin->contains( "pep_1" ) );
    REQUIRE( first_bin->contains( "pep_2" ) );
    REQUIRE( first_bin->contains( "pep_3" ) );
    REQUIRE( !( first_bin->contains( "pep_4" ) ) );

    auto second_bin = first_bin + 1;
    REQUIRE( second_bin->contains( "pep_4" ) );
    REQUIRE( second_bin->contains( "pep_5" ) );
    REQUIRE( second_bin->contains( "pep_6" ) );
    REQUIRE( !( second_bin->contains( "pep_9" ) ) );

	// test smallest() and add_bins() method
	first_bin->add_peptide("pep_10");
	first_bin->add_peptide("pep_11");
	first_bin->add_peptide("pep_12");

	REQUIRE(first_bin->size() == 6);
	REQUIRE(second_bin->size() == 3);
	REQUIRE(bin_c.smallest() == *second_bin);
}

TEST_CASE("", "[peptide_scoring]")
{
}

TEST_CASE( "Verify sums and minimum bin size resizing", "[module_bin]" )
{
    SECTION( "Summing Counts of items in a matrix" )
    {
        module_bin mod;
        std::vector<std::string> row_labels{ "row_1", "row_2", "row_3" };
        std::vector<std::string> col_labels{ "col_1", "col_2", "col_3" };

        labeled_matrix<double,std::string>
            matrix{ 3,
                    3,
                    row_labels,
                    col_labels
                };

        matrix.set_all( 4 );

        auto summed_counts = mod.sum_counts( matrix );

        REQUIRE( summed_counts[ 0 ] == 12 );
        REQUIRE( summed_counts[ 1 ] == 12 );
        REQUIRE( summed_counts[ 2 ] == 12 );

        matrix( 2, 2 ) = 11.5;
        summed_counts = mod.sum_counts( matrix );

        REQUIRE( summed_counts[ 2 ] == 19.5 );
        REQUIRE( summed_counts[ 0 ] == 12 );
    }

    SECTION( "Tweak peptide combination when a bin falls below minimum size" )
    {
        module_bin mod = module_bin();
        options_bin opt = options_bin();
        opt.min_bin_size = 300;
        opt.input_scores_fname = "../test/input_data/test_bins.tsv";
        opt.rounding_factor = 1;
        opt.output_bins_fname = "../test/test_bin_output.txt";
        mod.run( &opt );
        std::ifstream bin_output( opt.output_bins_fname, std::ios_base::in );
        std::ifstream bin_expected_output( "../test/expected/test_expected_bins.tsv",
                                                        std::ios_base::in );
        std::string expected_line;
        std::string output_line;

        while( std::getline( bin_output, output_line ) &&
                                    std::getline( bin_expected_output, expected_line ) )
            {
                std::unordered_set<std::string> bin_output_set;
                std::unordered_set<std::string> bin_expected_set;
                boost::split( bin_output_set, output_line, boost::is_any_of( "\t" ) );
                boost::split( bin_expected_set, expected_line, boost::is_any_of( "\t" ) );
                for( auto probe: bin_expected_set )
                    {
                        REQUIRE( bin_output_set.find( probe ) !=
                                                    bin_output_set.end() );
                    }
            }
    }
}

TEST_CASE( "Ranking Probes based upon their scores", "[probe_rank]" )
{
    std::unordered_map<std::string, double> probe_scores
    {
        { "p1", 1.445 },
        { "p2", 4.655 },
        { "p3", 1.045 },
        { "p4", 156.095 }
    };

    auto rank_probes = [&]( probe_rank& pr )
        {
            // we want to force an order here
            pr.rank_probe( 1.445, "p1" );
            pr.rank_probe( 4.655, "p2" );
            pr.rank_probe( 1.045, "p3" );
            pr.rank_probe( 156.095, "p4" );
        };

    SECTION( "Rounding to the nearest integer" )
        {
            probe_rank int_probe_rank{ 0 };
            rank_probes( int_probe_rank );
            decltype( int_probe_rank )::rank_track_type
                expected_values
            {
                { 1.0, { "p1", "p3" } },
                { 5.0, { "p2" } },
                { 156.0, { "p4" } },
            };

            REQUIRE( int_probe_rank.get_probe_ranks() == expected_values );
			REQUIRE(*int_probe_rank.get_probes_with_rank(1) == *expected_values.find(1.0));
			REQUIRE(*int_probe_rank.get_probes_with_rank(5.005) == *expected_values.find(5.0));
			REQUIRE(*int_probe_rank.get_probes_with_rank(156.056) == *expected_values.find(156.0));
        }
    SECTION( "Rounding to the tenth's place" )
        {
            probe_rank int_probe_rank{ 1 };
            rank_probes( int_probe_rank );
            decltype( int_probe_rank )::rank_track_type
                expected_values
            {
                { 1.0, { "p3" } },
                { 1.4, { "p1" } },
                { 4.7, { "p2" } },
                { 156.09999999999999, { "p4" } }, // hard-coded value, found by inspecting
                                                  // values in debugger.
            };

            REQUIRE( int_probe_rank.get_probe_ranks() == expected_values );
			REQUIRE(*int_probe_rank.get_probes_with_rank(4.7445) == *expected_values.find(4.7));
			REQUIRE(*int_probe_rank.get_probes_with_rank(1.4483) == *expected_values.find(1.4));
        }
    SECTION( "Rounding to the hundredth's place" )
        {
            probe_rank int_probe_rank{ 2 };
            rank_probes( int_probe_rank );
            decltype( int_probe_rank )::rank_track_type
                expected_values
            {
                { 1.05, { "p3" } },
                { 1.45, { "p1" } },
                { 4.66, { "p2" } },
                { 156.10, { "p4" } },
            };

            REQUIRE( int_probe_rank.get_probe_ranks() == expected_values );
			REQUIRE(*int_probe_rank.get_probes_with_rank(156.10) == *expected_values.find(156.10));
			REQUIRE(*int_probe_rank.get_probes_with_rank(1.05) == *expected_values.find(1.05));
			REQUIRE(*int_probe_rank.get_probes_with_rank(1.4522) == *expected_values.find(1.45));
        }
}

TEST_CASE( "Unary Predicate Reduction", "[module_enrich]" )
{
    using namespace predicate;

    auto positive = []( const int x ) -> bool { return x > 0; };
    auto gt_10 = []( const int x ) -> bool { return x > 10; };
    auto lt_100 = []( const int x ) -> bool { return x < 100; };
    auto odd = []( const int x ) -> bool { return x % 2; };

    // assert 5 is positive
    REQUIRE( value_constrained_by( 5, positive ) );

    // -5 is not positive
    REQUIRE( !value_constrained_by( -5, positive ) );

    // 51 is positive, greater than 10, less than 100,
    // and odd
    REQUIRE( value_constrained_by( 51,
                                   positive,
                                   gt_10,
                                   lt_100,
                                   odd
                                  )
             );

    // 51 is not odd
    REQUIRE( !value_constrained_by( 50,
                                    positive,
                                    gt_10,
                                    lt_100,
                                    odd
                                  )
             );

}

TEST_CASE( "Valid For", "[module_enrich]" )
{
    std::vector<int> vals{ 1, 2, 3, 4, 5 };

    auto lt_100 = []( const int x ) -> bool { return x < 100; };
    auto positive = []( const int x ) -> bool { return x > 0; };
    auto even = []( const int x ) -> bool { return x % 2 == 0; };
    using namespace predicate;

    SECTION( "Default function, no use of special 'get' function" )

        {
            auto dest = valid_for( vals.begin(),
                                   vals.end(),
                                   vals.begin(),
                                   lt_100
                                 );

            REQUIRE( dest == vals.end() );

            std::vector<int> vals2{ 1, -1, 3, 4, 5 };
            std::vector<int> vals_result{};
            valid_for( vals2.begin(),
                       vals2.end(),
                       std::back_inserter( vals_result ),
                       positive
                       );

            // the expression will have been false for one value
            REQUIRE( std::distance( vals_result.begin(), vals_result.end() ) == 4 );
            REQUIRE( std::find( vals_result.begin(), vals_result.end(), -1 ) == vals_result.end() );
        }

    SECTION( "Use of custom get function" )
        {
            using pair_t = std::pair<int,std::string>;
            std::vector<pair_t>
                values{
                        { 5, "str1" },
                        { 6, "str2" },
                        { 9, "str3" },
                        { 1, "str4" },
                       };
            auto get_first = []( const pair_t& val )
                -> double
                {
                    return val.first;
                };

            std::vector<pair_t> result;
            std::vector<pair_t> expected{ { 6, "str2" } };

            valid_for( values.begin(),
                       values.end(),
                       std::back_inserter( result ),
                       even,
                       get_first
                     );

            REQUIRE( result.size() == 1 );
            REQUIRE( result == expected );
        }
}

TEST_CASE( "Use of predicate logic", "[predicate]" )
{
    bool p = true;
    bool q = false;

    REQUIRE( !predicate::implication( p, q ) );
    REQUIRE(  predicate::implication( q, p ) );

    REQUIRE( predicate::disjunction( p, q ) );
    REQUIRE( predicate::disjunction( q, p ) );

    REQUIRE( !predicate::conjunction( p, q ) );
    REQUIRE( !predicate::conjunction( q, p ) );

    REQUIRE( predicate::biconditional( !p, q ) );
    REQUIRE( predicate::biconditional( !q, p ) );

}

TEST_CASE( "Parsing samples file", "[module_enrich]" )
{
    module_enrich mod;
    std::string input_str;
    input_str = "sample1\tsample2\nsample3\tsample4\n";

    std::istringstream file;
    file.str( input_str );
    auto samples = mod.parse_samples( file );

    REQUIRE( samples.size() == 2 );

    file.clear();

    file.str( "sample1\tsample2\tsample3\nsample4\tsample5\nsample6\t\tsample7" );

    samples = mod.parse_samples( file );

    REQUIRE( samples.size() == 3 );
}

TEST_CASE( "Meeting the threshold for a pair", "[module_enrich]" )
{
    module_enrich mod;
    auto mp = []( const int a, const int b ) -> std::pair<int,int>
        {
            return std::make_pair( a, b );
        };

    bool met = mod.pair_threshold_met(  mp( 1, 2 ),
                                        mp( 3, 4 )
                                     );
    REQUIRE( !met );

    met = mod.pair_threshold_met(  mp( 1, 2 ),
                                   mp( 2, 1 )
                                );
    REQUIRE( met );

    met = mod.pair_threshold_met(  mp( 15, 4 ),
                                   mp( 16, 3 )
                                );
    REQUIRE( !met );

    met = mod.pair_threshold_met(  mp( 10, 20 ),
                                   mp( 3, 4 )
                                );
    REQUIRE( met );

}

TEST_CASE("Test enrich drops replicates with low scores if --low_raw_reads flag passed", "[module_enrich]")
{
	// initialize enrich components
	module_enrich e_mod;
	options_enrich e_opts;

	std::pair<std::string, std::string> pair1{"../test/input_data/test_enrich_Z-HDI75.tsv", "10"};
	std::pair<std::string, std::string> pair2{"../test/input_data/test_enrich_CS.tsv", "20"};

	e_opts.matrix_thresh_fname_pairs = {pair1, pair2};
	e_opts.in_samples_fname = "../test/input_data/test_enrich_PN.tsv";
	e_opts.in_raw_scores_fname = "../test/input_data/test_enrich_raw_scores.tsv";
	e_opts.raw_scores_params_str = "70";
	e_opts.out_suffix = "_enriched_output.txt";
	e_opts.out_enrichment_failure = "test_enrich_fail_output.txt";
	e_opts.out_dirname = "../test/test_enrich_Z-HDI75_70raw_output";
	e_opts.low_raw_reads = true;

	// run enrich
	e_mod.run(&e_opts);
	
	// check failed enrichment file for properly formatted log of dropped replicates
	std::ifstream ifexpected(
		"../test/expected/test_expected_enrich_Z-HDI75_70raw/test_expected_enrich_fail.tsv",
		std::ios_base::in
	);
	std::ifstream ifactual(
		"../test/test_enrich_Z-HDI75_70raw_output/test_enrich_fail_output.txt",
		std::ios_base::in
	);
	std::string expected_line = "";
	std::string actual_line = "";

	while (!ifexpected.eof() && !ifactual.eof())
	{
		std::getline(ifexpected, expected_line);
		std::getline(ifactual, actual_line);
		REQUIRE(expected_line.compare(actual_line) == 0);
	}
}

TEST_CASE( "File IO read_file function", "[file_io]" )
{
    std::istringstream stream;
    stream.str( "First\t0.45\nSecond\t0.58\n" );
    auto create_fn = []( typename std::vector<std::string>::iterator begin,
                         typename std::vector<std::string>::iterator end
                         ) -> std::pair<std::string,double>
        {
            std::string first = *begin;
            std::string second = *(++begin);
            double second_double = 0.0;

            if( begin != end )
                {
                    second_double = std::strtod( second.c_str(), nullptr );
                }

            return std::make_pair( first, second_double );
        };

    std::vector<std::pair<std::string,double>> values;

    pepsirf_io::read_file( stream,
                           boost::is_any_of( "\t" ),
                           create_fn,
                           std::back_inserter( values )
                         );

    REQUIRE( values.size() == 2 );
    REQUIRE( values[ 0 ] == std::pair<std::string,double>( "First", 0.45 ) );
    REQUIRE( values[ 1 ] == std::pair<std::string,double>( "Second", 0.58 ) );

    stream.str( "" );
    stream.clear();
    values.clear();

    stream.str( "First\t0.45\nSecond\t0.58" );

    pepsirf_io::read_file( stream,
                           boost::is_any_of( "\t" ),
                           create_fn,
                           std::back_inserter( values )
                         );

    REQUIRE( values.size() == 2 );
    REQUIRE( values[ 0 ] == std::pair<std::string,double>( "First", 0.45 ) );
    REQUIRE( values[ 1 ] == std::pair<std::string,double>( "Second", 0.58 ) );


}

TEST_CASE( "Testing Codon -> AA maps", "[nt_aa_map]" )
{

    auto mapped = codon_aa_mappings::default_codon_aa_map( "TCC" );
    REQUIRE( mapped == 'S' );

    mapped = codon_aa_mappings::default_codon_aa_map( "GTC" );

    REQUIRE( mapped == 'V' );
}

TEST_CASE( "Testing nt->aa translation", "[nt_aa_translator]" )
{
    nt_aa_translator<codon_aa_mappings::default_map_type>
        translator( codon_aa_mappings::default_codon_aa_map );

    sequence nt_seq{ "Seq1", "AAAATACCC" };

    sequence translated = translator( nt_seq );

    REQUIRE( translated.seq == "KIP" );

    nt_seq.seq = "TAGTAATGATTT";
    translated = translator( nt_seq );

    REQUIRE( translated.seq == "___F" );

}

#ifdef ZLIB_ENABLED

TEST_CASE( "Reading/Writing Gzipped information", "[pepsirf_io]" )
{
    std::string data = "One banana, two banana, three banana, split!";


    // I'm not certain why this needs its own block.
    // Trying to close/flush the output stream manually does not work,
    // but allowing it to be closed by the destructor DOES
    // This is only a problem when reading/writing in the same
    // breath, typically only reading or writing will be done, and this will
    // not be a problem.
    {
        std::ofstream out{ "test.gz", std::ios_base::out | std::ios_base::binary };
        pepsirf_io::gzip_writer write{ out };
        write << data;
    }

    std::ifstream in{ "test.gz",
                      std::ios_base::in | std::ios_base::binary
                    };

    pepsirf_io::gzip_reader read{ in };

    std::string comp;

    std::getline( read, comp );

    REQUIRE( comp == data );
}

#endif // ZLIB_ENABLED

TEST_CASE( "Determining whether a file is gzipped.", "[pepsirf_io]" )
{
    std::ifstream true_expected{ "test.gz" };

    REQUIRE( pepsirf_io::is_gzipped( true_expected ) );

    std::ifstream false_expected{ "../test/input_data/test.fasta" };
    REQUIRE( !pepsirf_io::is_gzipped( false_expected ) );
}

TEST_CASE( "Subjoin name list filter is optional", "[module_subjoin]" )
{
    module_subjoin mod = module_subjoin();
    options_subjoin opts = options_subjoin();
    opts.use_sample_names = true;
    opts.out_matrix_fname = "../test/test_subjoin_output.txt";
    opts.input_matrix_name_pairs.emplace_back( std::make_pair( "../test/input_data/test_score_matrix.tsv", "" ) );
    mod.run( &opts );
}

TEST_CASE("Full test of link module's individual funcitons", "module_link")
{
	SECTION("Testing prot map creation")
	{
	}
	SECTION("Testing pep map creation")
	{
	}
	SECTION("Testing pep map creation with kmer penalty")
	{
	}
}

TEST_CASE( "Metadata file can be given in place of taxonomic id index", "[module_link]" )
{
    SECTION( "Verifying map creation and target sequence retrieval individually from link module" )
    {
        std::ifstream metadata_file( "../test/test_meta_fraction.metadata", std::ios_base::in );
        module_link mod = module_link();
        std::unordered_map<std::string, std::string> meta_map(
                                                    { { ">ID=A0A2Z4GZU5_HHV1 AC=A0A2Z4GZU5 OXX=10298,10298,10294,10292", "10294" },
                                                    { ">ID=A0A2L2Q9I3_9BETA AC=A0A2L2Q9I3 OXX=32603,32603,40272,10292", "40272" },
                                                    { ">ID=E9RHT2_HPV67 AC=E9RHT2 OXX=37120,337041,333750,151340", "333750" } }
            );
        int id_index = 2;
        std::string id_index_spec_id = mod.get_id( ">ID=A0A2Z4GZU5_HHV1 AC=A0A2Z4GZU5 OXX=10298,10298,10294,10292", id_index );
        std::string metadata_spec_id = metadata_map::get_id( ">ID=A0A2Z4GZU5_HHV1 AC=A0A2Z4GZU5 OXX=10298,10298,10294,10292", &meta_map );
        REQUIRE( id_index_spec_id.compare( metadata_spec_id ) == 0 );
    }
    SECTION( "Verifying metadata file usage feature successfully performs and integrates with link module." )
    {
        module_link mod;
        options_link opts;
        opts.metadata_fname =  "../test/input_data/full_design_clean_min30_taxtweak_100perc_jingmens_2019-09-12_segment.metadata,Name,Species";
        opts.prot_file_fname = "../test/input_data/full_design_clean_min30_taxtweak_100perc_jingmens_2019-09-12_segment.fasta";
        opts.peptide_file_fname = "../test/input_data/PV1_10K3000_53_encoded_segment.faa";
        opts.kmer_size = 7;
        opts.output_fname = "../test/test_metadata_output.txt";
        mod.run( &opts );
        REQUIRE( !opts.output_fname.empty() );
    }
}

TEST_CASE( "Discard outliers of bin distribution by trimming for further calculation", "[module_zscore]" )
{
    SECTION( "Verify highest density interval option properly trims distribution" )
    {
        module_zscore mod = module_zscore();
        options_zscore opts = options_zscore();
        peptide_score_data_sample_major input;
        opts.in_fname = "../../../../Downloads/combo_norm_avgSBdiff.tsv";
        opts.in_bins_fname = "../../../../Downloads/NS30_IgG_combo_SB4_300r1.tsv";
        opts.out_fname = "../../../../Downloads/zscore_test_hdi_output.tsv";
        opts.hpd_percent = 0.75;
        peptide_scoring::parse_peptide_scores( input, opts.in_fname );
        std::ifstream bins_file( opts.in_bins_fname, std::ios_base::in );
        bin_collection peptide_bins = peptide_bin_io::parse_bins( bins_file );
        auto& zscore_matrix = input.scores;
        std::vector<nan_report> nan_values;
        for( const auto sample_pair : zscore_matrix.get_row_labels() )
            {
                std::string sample_name = sample_pair.first;

                mod.calculate_zscores(  peptide_bins,
                                        &opts,
                                        zscore_matrix,
                                        sample_name,
                                        std::back_inserter( nan_values )
                                     );
            }
        peptide_scoring::write_peptide_scores( opts.out_fname, input );
    }
    SECTION( "Verify symmetrical trim option properly trims both head and tail" )
    {
        module_zscore mod = module_zscore();
        options_zscore opts = options_zscore();
        peptide_score_data_sample_major input;
        opts.in_fname = "../../../../Downloads/combo_norm_avgSBdiff.tsv";
        opts.in_bins_fname = "../../../../Downloads/NS30_IgG_combo_SB4_300r1.tsv";
        opts.out_fname = "../../../../Downloads/zscore_test_trim_output.tsv";
        opts.hpd_percent = 0.0;
        opts.trim_percent = 0.90;
        peptide_scoring::parse_peptide_scores( input, opts.in_fname );
        std::ifstream bins_file( opts.in_bins_fname, std::ios_base::in );
        bin_collection peptide_bins = peptide_bin_io::parse_bins( bins_file );
        auto& zscore_matrix = input.scores;
        std::vector<nan_report> nan_values;
        for( const auto sample_pair : zscore_matrix.get_row_labels() )
            {
                std::string sample_name = sample_pair.first;

                mod.calculate_zscores(  peptide_bins,
                                        &opts,
                                        zscore_matrix,
                                        sample_name,
                                        std::back_inserter( nan_values )
                                     );
            }
        peptide_scoring::write_peptide_scores( opts.out_fname, input );
    }
}
