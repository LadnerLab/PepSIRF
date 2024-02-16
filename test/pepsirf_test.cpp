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

#include "bk_tree.h"
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

	sequence seq2("sequence2", "ATGCGGGGTC");

	std::string new_seq_name = "sequence1";
	seq.set_name(new_seq_name);
	REQUIRE(seq.name.compare("sequence1") == 0);
	REQUIRE(seq == seq2);

	REQUIRE(seq.length() == 10);
	REQUIRE(seq2.length() == seq.length());
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
	SECTION("set_intersection() with Get")
	{
		auto only_positives = [](int num)
		{
			return num > 0 ? num : 0;
		};

		std::set<int> int_set1 = {10, -13, 3, -4, -5, 4};
		std::set<int> int_set2 = {10, -1, 3, 4, 5, 8};
		std::set<int> dest_set;
		std::vector<int> expected_set = {
			3, 4, 10
		};

		set_intersection(dest_set, int_set1, int_set2, only_positives);

		REQUIRE(dest_set.size() == expected_set.size());
		std::size_t i = 0;
		for (auto dest_it : dest_set)
		{
			REQUIRE(dest_it == expected_set[i]);
			i += 1;
		}
	}
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
	SECTION("set_union() with Get")
	{
		std::unordered_map<char, int> map1 = {
			{'c', 8}, {'b', 10}, {'y', 3}, {'w', 1}
		};
		std::unordered_map<char, int> map2 = {
			{'m', 13}, {'q', 12}, {'o', 0}, {'z', 3}
		};
		std::set<int> dest_set;
		std::vector<int> expected_set = {
			'b', 'c', 'm', 'o', 'q', 'w', 'y', 'z'
		};

		set_union(dest_set, map1, map2, get_key<char, int>());

		REQUIRE(dest_set.size() == expected_set.size());
		std::size_t i = 0;
		for (auto dest_it : dest_set)
		{
			REQUIRE(dest_it == expected_set[i]);
			i += 1;
		}
	}
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

        std::string expected_line = "";
        std::string actual_line = "";
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

        std::string expected_line;
        std::string actual_line;
        while (!ifexpected.eof() && !ifactual.eof())
        {
            std::getline(ifexpected, expected_line);
            std::getline(ifactual, actual_line);
            REQUIRE(expected_line.compare(actual_line) == 0);
        }
	}
    SECTION("info module properly calculates and formats column sums file")
    {
        std::string line;
        std::set<std::string> expected_set;
        std::set<std::string> actual_set;

        // collect expected lines into a set
        {
            std::ifstream ifexpected(
                "../test/expected/test_expected_info_col_sums.tsv",
                std::ios_base::in
            );

            while (!ifexpected.eof())
            {
                std::getline(ifexpected, line);
                expected_set.insert(line);
            }
        }
        // collect actual lines into a set
        {
            std::ifstream ifactual(
                "../test/test_info_col_sums.tsv",
                std::ios_base::in
            );

            while (!ifactual.eof())
            {
                std::getline(ifactual, line);
                actual_set.insert(line);
            }
        }

        // compare expected and actual sets
        REQUIRE(actual_set.size() == expected_set.size());

        auto expected_ref = expected_set.find("BB.5_1X_NS30_B\t2995.00");
        REQUIRE(actual_set.find(*expected_ref) != actual_set.end());
        
        expected_ref = expected_set.find("BB.6_1X_NS30_A\t49.00");
        REQUIRE(actual_set.find(*expected_ref) != actual_set.end());

        expected_ref = expected_set.find("BB.3_1X_NS30_B\t1604.00");
        REQUIRE(actual_set.find(*expected_ref) != actual_set.end());

        expected_ref = expected_set.find("BB.5_1X_NS30_C\t1801.00");
        REQUIRE(actual_set.find(*expected_ref) != actual_set.end());

        expected_ref = expected_set.find("BB.1_1X_NS30_A\t327.00");
        REQUIRE(actual_set.find(*expected_ref) != actual_set.end());

        expected_ref = expected_set.find("BB.3_1X_NS30_A\t1251.00");
        REQUIRE(actual_set.find(*expected_ref) != actual_set.end());

        expected_ref = expected_set.find("BB.5_1X_NS30_A\t1597.00");
        REQUIRE(actual_set.find(*expected_ref) != actual_set.end());

        expected_ref = expected_set.find("BB.6_1X_NS30_B\t963.00");
        REQUIRE(actual_set.find(*expected_ref) != actual_set.end());
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

        std::string expected_line;
        std::string actual_line;

        while (!ifexpected.eof() && !ifactual.eof())
        {
            std::getline(ifexpected, expected_line);
            std::getline(ifactual, actual_line);
            REQUIRE(expected_line.compare(actual_line) == 0);
        }
	}
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

TEST_CASE("Full test of util's individual methods", "[util]")
{
	SECTION("Test is_integer()")
	{
		SECTION("Is 9 an integer")
		{
			REQUIRE(is_integer(9));
		}
		SECTION("Is 6.83 an integer")
		{
			REQUIRE(!is_integer(6.38));
		}
	}

	SECTION("Test difference operation")
	{
		SECTION("Difference between 22 and 11")
		{
			REQUIRE(difference<int>()(22, 11) == 11);
		}
		SECTION("Difference between 2546.581 and 64.32648")
		{
			REQUIRE(difference<double>()(2546.581, 64.32648) == 2482.25452);
		}
		SECTION("Difference between 100.00 and 1000.00")
		{
			REQUIRE(difference<double>()(100.00, 1000.00) == -900.00);
		}
	}

	SECTION("Test ratio operation")
	{
		SECTION("Ratio of 1 and 3")
		{
			REQUIRE(ratio<double>()(1, 3) == 1.0 / 3.0);
		}
		SECTION("Ratio of 10 and 20")
		{
			REQUIRE(ratio<double>()(20, 10) == 1.0 / 2.0);
		}
		SECTION("Ration of 100 and 38")
		{
			REQUIRE(ratio<double>()(100, 38) == 38.0 / 100.0);
		}
	}

	SECTION("Test pair positional comparison")
	{
		auto both_even = [](int num, int other)
		{
			return (num % 2 == 0) && (other % 2 == 0);
		};
		auto both_odd = [](int num, int other)
		{
			return (num % 2 == 1) && (other % 2 == 1);
		};
		auto both_positive = [](int num, int other)
		{
			return (num > 0) && (other > 0);
		};

		std::pair<int, int> first = {10, 32};
		std::pair<int, int> second = {6, 18};

		REQUIRE(pair_positional_compare(
				first, second,
				both_even, both_odd
			) == false
		);

		first = {8, 13};
		second = {22, 583};

		REQUIRE(pair_positional_compare(
				first, second,
				both_even, both_odd
			) == true
		);

		first = {3, 11};
		second = {21, 5};

		REQUIRE(pair_positional_compare(
				first, second,
				both_odd, both_odd
			) == true
		);

		first = {-10, 44};
		second = {12, 0};

		REQUIRE(pair_positional_compare(
				first, second,
				both_positive, both_positive
			) == false
		);
	}

	SECTION("Test non-increasing pair comparison")
	{
		std::vector<std::pair<char, char>> pair_vec = {
			{'c', 'a'}, {'z', 's'}, {'n', 'z'}, {'w', 't'}, {'b', 'o'}
		};
		std::vector<std::pair<char, char>> expected_pair_vec = {
			{'n', 'z'}, {'w', 't'}, {'z', 's'}, {'b', 'o'}, {'c', 'a'}
		};

		std::sort(
			pair_vec.begin(), pair_vec.end(),
			compare_pair_non_increasing<char, char>()
		);

		for (std::size_t i = 0; i < expected_pair_vec.size(); i += 1)
		{
			REQUIRE(pair_vec[i].second == expected_pair_vec[i].second);
		}
	}

	SECTION("Test non-decreasing pair comparison")
	{
		std::vector<std::pair<char, char>> pair_vec = {
			{'c', 'a'}, {'z', 's'}, {'n', 'z'}, {'w', 't'}, {'b', 'o'}
		};
		std::vector<std::pair<char, char>> expected_pair_vec = {
			{'c', 'a'}, {'b', 'o'}, {'z', 's'}, {'w', 't'}, {'n', 'z'}
		};

		std::sort(
			pair_vec.begin(), pair_vec.end(),
			compare_pair_non_decreasing<char, char>()
		);

		for (std::size_t i = 0; i < expected_pair_vec.size(); i += 1)
		{
			REQUIRE(pair_vec[i].second == expected_pair_vec[i].second);
		}
	}

	SECTION("Test conversion of pairs to map")
	{
		std::vector<std::pair<std::string, int>> my_pair_vec = {
			{"pep1", 10}, {"pep2", 22}, {"pep52", 0}, {"pep100", 120}
		};
		std::map<std::string, int> my_map_of_pairs;

		pairs_to_map(my_map_of_pairs, my_pair_vec.begin(), my_pair_vec.end());

		auto it = my_map_of_pairs.find(my_pair_vec[0].first);
		REQUIRE(it->second == 10);
		it = my_map_of_pairs.find(my_pair_vec[1].first);
		REQUIRE(it->second == 22);
		it = my_map_of_pairs.find(my_pair_vec[2].first);
		REQUIRE(it->second == 0);
		it = my_map_of_pairs.find(my_pair_vec[3].first);
		REQUIRE(it->second == 120);
	}

	SECTION("Test container printing")
	{
		std::ostringstream container_print_out;
		std::string expected_print_out = "[pep1, 7]\n[pep2, 3]\n[pep3, 10]\n[pep4, 6]\n";
		std::vector<std::pair<std::string, int>> container = {
			{"pep1", 7}, {"pep2", 3}, {"pep3", 10}, {"pep4", 6}
		};

		print_container(container_print_out, container, "\n");

		REQUIRE(expected_print_out.compare(container_print_out.str()) == 0);
	}
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

TEST_CASE("pairwise distance symmetry optimized", "[distance_tools]")
{
    distance_matrix<int> dist(5);
	std::vector<int> int_vec = {6, 8, 2, 0, 1};

	auto dist_calc = [](const int a, const int b)
	{
		return a / std::abs(a - b);
	};

	REQUIRE(dist_calc(4, 8) == 1);

	distance_tools::pairwise_distances_symmetry_optimized(
		dist, int_vec.begin(), int_vec.end(), dist_calc
	);
	
	REQUIRE(dist[0][0] == 3);
	REQUIRE(dist[0][1] == 1);
	REQUIRE(dist[0][2] == 1);
	REQUIRE(dist[0][3] == 1);
	REQUIRE(dist[1][0] == 1);
	REQUIRE(dist[1][1] == 1);
	REQUIRE(dist[1][2] == 1);
	REQUIRE(dist[2][0] == 1);
	REQUIRE(dist[2][1] == 2);
	REQUIRE(dist[3][0] == 0);
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

	SECTION("deconv orders tied species by ID in output file")
	{
		std::ifstream ifexpected(
			"../test/expected/test_expected_deconv_output.tsv",
			std::ios_base::in
		);
		std::ifstream ifactual(
			"../test/test_deconv_output.tsv",
			std::ios_base::in
		);
		std::string expected_line = "";
		std::string actual_line = "";

		while(!ifexpected.eof() && !ifactual.eof())
			{
				std::getline(ifexpected, expected_line);
				std::getline(ifactual, actual_line);
				REQUIRE(expected_line.compare(actual_line) == 0);
			}
	}
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

	REQUIRE(pep.get_name() == "pep1");
	REQUIRE(pep2.get_name() == "p11");

    REQUIRE( pep == pep2 );

    pep2.set_sequence( "AGGG" );

    REQUIRE( pep != pep2 );
    REQUIRE( !pep2.get_sequence().compare( "AGGG" ) );

	pep.set_name("p11");

	REQUIRE(pep.get_name() == pep2.get_name());
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

     std::unordered_map<std::string, scored_peptide<double>> highest_scores;

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
		std::vector<double> dest_vec = {0, 0, 0};
		matrix<double> doubles_mat(3, 3);

		doubles_mat(0, 0) = -9.462;
		doubles_mat(1, 0) = -10.462;
		doubles_mat(2, 0) = -0.492;
		
		doubles_mat(0, 1) = 5.632;
		doubles_mat(1, 1) = -0.843;
		doubles_mat(2, 1) = 0.748;
		
		doubles_mat(0, 2) = 2.351;
		doubles_mat(1, 2) = 3.245;
		doubles_mat(2, 2) = 0.134;

		norm_mod.get_sum(&dest_vec, &doubles_mat);

		REQUIRE(dest_vec[0] == -20.416);
		REQUIRE(dest_vec[1] == 5.537);
		REQUIRE(dest_vec[2] == 5.730);
	}
	SECTION("Test get_neg_average() works as expected")
	{
		peptide_score_data score_data;
		score_data.pep_names = {"pep1", "pep2", "pep3", "pep4"};
		score_data.sample_names = {"sample1", "sample2", "sample3", "sample4"};

		score_data.scores = labeled_matrix<double, std::string>(
			score_data.pep_names.size(), score_data.sample_names.size()
		);

		score_data.scores(0, 0) = 6.388;
		score_data.scores(0, 1) = 10.997;
		score_data.scores(0, 2) = 0.037;
		score_data.scores(0, 3) = 5.730;

		score_data.scores(1, 0) = 38.133;
		score_data.scores(1, 1) = -7.483;
		score_data.scores(1, 2) = 4.909;
		score_data.scores(1, 3) = 12.283;
		
		score_data.scores(2, 0) = -3.182;
		score_data.scores(2, 1) = 0.291;
		score_data.scores(2, 2) = -0.321;
		score_data.scores(2, 3) = -20.441;

		score_data.scores(3, 0) = -0.210;
		score_data.scores(3, 1) = 4.008;
		score_data.scores(3, 2) = -0.110;
		score_data.scores(3, 3) = -1.001;

		std::unordered_set<std::string> neg_filter = {"sample3", "sample4"};
		std::unordered_map<std::string, double> dest_pep_avgs;

		norm_mod.get_neg_average(&score_data, &neg_filter, &dest_pep_avgs);

		auto pep_avg = dest_pep_avgs.find("pep1");
		REQUIRE(pep_avg->second == 2.8835);

		pep_avg = dest_pep_avgs.find("pep2");
		REQUIRE(pep_avg->second == 8.596);

		pep_avg = dest_pep_avgs.find("pep3");
		REQUIRE(pep_avg->second == -10.381);

		pep_avg = dest_pep_avgs.find("pep4");
		REQUIRE(pep_avg->second == -0.5555);
	}
	SECTION("Test constant_factor_normalization() works as expected")
	{
		SECTION("Normalize column sums by 2")
		{
			std::vector<double> columns = {
				6.832, -0.003, 10.374, 12.889, 1.234
			};

			norm_mod.constant_factor_normalization(&columns, 2);

			REQUIRE(columns[0] == 3.416);
			REQUIRE(columns[1] == -0.0015);
			REQUIRE(columns[2] == 5.187);
			REQUIRE(columns[3] == 6.4445);
			REQUIRE(columns[4] == 0.617);
		}
		SECTION("Normalize column sums by 4")
		{
			std::vector<double> columns = {
				6.832, -0.003, 10.374, 12.889, 1.234
			};

			norm_mod.constant_factor_normalization(&columns, 4);

			REQUIRE(columns[0] == 1.708);
			REQUIRE(columns[1] == -0.00075);
			REQUIRE(columns[2] == 2.5935);
			REQUIRE(columns[3] == 3.22225);
			REQUIRE(columns[4] == 0.3085);
		}
		SECTION("Normalize column sums by 100")
		{
			std::vector<double> columns = {
				16.832, -15.033, 10.374, 12.889, 10.234
			};

			norm_mod.constant_factor_normalization(&columns, 100);

			REQUIRE(columns[0] == 0.16832);
			REQUIRE(columns[1] == -0.15033);
			REQUIRE(columns[2] == 0.10374);
			REQUIRE(columns[3] == 0.12889);
			REQUIRE(columns[4] == 0.10234);
		}
	}
	SECTION("Test filter_neg_control_start() works as expected")
	{
		peptide_score_data score_data;
		score_data.sample_names = {"GYT_UU_C22", "NAE_COV_B45", "NAE_COV_A12", "CNB_MO8"};

		std::unordered_set<std::string> neg_filter;

		norm_mod.filter_neg_control_start(&score_data, &neg_filter, "NAE");

		REQUIRE(neg_filter.find("NAE_COV_B45") != neg_filter.end());
		REQUIRE(neg_filter.find("NAE_COV_A12") != neg_filter.end());
	}
	SECTION("Test normalize_counts() works as expected")
	{
		matrix<double> score_mat(4, 4);

		score_mat(0, 0) = 5.683;
		score_mat(0, 1) = -10.483;
		score_mat(0, 2) = 3.000;

		score_mat(1, 0) = 12.489;
		score_mat(1, 1) = 39.542;
		score_mat(1, 2) = 10.000;

		score_mat(2, 0) = -0.004;
		score_mat(2, 1) = -1.834;
		score_mat(2, 2) = 0.100;

		score_mat(3, 0) = 0.032;
		score_mat(3, 1) = -1.002;
		score_mat(3, 2) = 6.000;
		
		std::vector<double> norm_factors = {
			1, 5, 0.5
		};

		norm_mod.normalize_counts(&score_mat, &norm_factors);

		REQUIRE(score_mat(0, 0) == 5.683);
		REQUIRE(score_mat(0, 1) == -2.0966);
		REQUIRE(score_mat(0, 2) == 6.000);

		REQUIRE(score_mat(1, 0) == 12.489);
		REQUIRE(score_mat(1, 1) == 7.9084);
		REQUIRE(score_mat(1, 2) == 20.000);

		REQUIRE(score_mat(2, 0) == -0.004);
		REQUIRE(score_mat(2, 1) == -0.3668);
		REQUIRE(score_mat(2, 2) == 0.200);

		REQUIRE(score_mat(3, 0) == 0.032);
		REQUIRE(score_mat(3, 1) == -0.2004);
		REQUIRE(score_mat(3, 2) == 12.000);
	}
	SECTION("Test compute_size_factors() works as expected")
	{
		std::vector<double> dest_vec;
		matrix<double> counts(3, 3);

		counts(0, 0) = 700.000;
		counts(0, 1) = 189.000;
		counts(0, 2) = 62.000;

		counts(1, 0) = 5.000;
		counts(1, 1) = 4.000;
		counts(1, 2) = 10.000;

		counts(2, 0) = 3.000;
		counts(2, 1) = 2.000;
		counts(2, 2) = 0.111;

		norm_mod.compute_size_factors(&dest_vec, &counts);

		// compares out to nine decimal places without having to worry about
		// rounding errors
		REQUIRE((int)(dest_vec[0] * 1000000) == 3435288);
		REQUIRE((int)(dest_vec[1] * 1000000) == 937154);
		REQUIRE((int)(dest_vec[2] * 1000000) == 307426);
	}
	SECTION("Test compute_diff() works as expected")
	{
		peptide_score_data_sample_major norm_score_diffs;
		norm_score_diffs.sample_names = {"sample1", "sample2", "sample3"};
		norm_score_diffs.pep_names = {"pep1", "pep2", "pep3"};

		norm_score_diffs.scores = labeled_matrix<double, std::string>(
			norm_score_diffs.pep_names.size(),
			norm_score_diffs.sample_names.size(),
			norm_score_diffs.pep_names, norm_score_diffs.sample_names
		);

		norm_score_diffs.scores(0, 0) = 7.000;
		norm_score_diffs.scores(0, 1) = -1.000;
		norm_score_diffs.scores(0, 2) = 12.000;

		norm_score_diffs.scores(1, 0) = 9.000;
		norm_score_diffs.scores(1, 1) = 1.000;
		norm_score_diffs.scores(1, 2) = -8.000;
		
		norm_score_diffs.scores(2, 0) = 10.000;
		norm_score_diffs.scores(2, 1) = -5.000;
		norm_score_diffs.scores(2, 2) = 4.000;

		std::unordered_map<std::string, double> pep_avgs = {
			{"pep1", 10.000}, {"pep2", 0.300}, {"pep3", 5.000}
		};

		norm_mod.compute_diff(&norm_score_diffs, &pep_avgs);

		REQUIRE(norm_score_diffs.scores(0, 0) == -3.000);
		REQUIRE(norm_score_diffs.scores(0, 1) == -11.000);
		REQUIRE(norm_score_diffs.scores(0, 2) == 2.000);
		
		REQUIRE(norm_score_diffs.scores(1, 0) == 8.700);
		REQUIRE(norm_score_diffs.scores(1, 1) == 0.700);
		REQUIRE(norm_score_diffs.scores(1, 2) == -8.300);

		REQUIRE(norm_score_diffs.scores(2, 0) == 5.000);
		REQUIRE(norm_score_diffs.scores(2, 1) == -10.000);
		REQUIRE(norm_score_diffs.scores(2, 2) == -1.000);
	}
	SECTION("Test compute_diff_ratio() works as expected")
	{
		peptide_score_data norm_score_diff_ratios;
		norm_score_diff_ratios.sample_names = {
			"sample1", "sample2", "sample3"
		};
		norm_score_diff_ratios.pep_names = {"pep1", "pep2", "pep3"};

		norm_score_diff_ratios.scores = labeled_matrix<double, std::string>(
			norm_score_diff_ratios.pep_names.size(),
			norm_score_diff_ratios.sample_names.size(),
			norm_score_diff_ratios.pep_names,
			norm_score_diff_ratios.sample_names
		);

		norm_score_diff_ratios.scores(0, 0) = 4.000;
		norm_score_diff_ratios.scores(0, 1) = -2.000;
		norm_score_diff_ratios.scores(0, 2) = 12.000;

		norm_score_diff_ratios.scores(1, 0) = 9.000;
		norm_score_diff_ratios.scores(1, 1) = 1.000;
		norm_score_diff_ratios.scores(1, 2) = -8.000;
		
		norm_score_diff_ratios.scores(2, 0) = 10.000;
		norm_score_diff_ratios.scores(2, 1) = -5.000;
		norm_score_diff_ratios.scores(2, 2) = 4.000;

		std::unordered_map<std::string, double> pep_avgs = {
			{"pep1", 2.000}, {"pep2", 3.000}, {"pep3", -8.000}
		};

		norm_mod.compute_diff_ratio(&norm_score_diff_ratios, &pep_avgs);

		REQUIRE(norm_score_diff_ratios.scores("pep1", "sample1") == 1.000);
		REQUIRE(norm_score_diff_ratios.scores("pep1", "sample2") == -2.000);
		REQUIRE(norm_score_diff_ratios.scores("pep1", "sample3") == 5.000);

		REQUIRE(norm_score_diff_ratios.scores("pep2", "sample1") == 2.000);
		REQUIRE(norm_score_diff_ratios.scores("pep2", "sample2") == (-2.000 / 3.000));
		REQUIRE(norm_score_diff_ratios.scores("pep2", "sample3") == (-11.000 / 3.000));

		REQUIRE(norm_score_diff_ratios.scores("pep3", "sample1") == -2.250);
		REQUIRE(norm_score_diff_ratios.scores("pep3", "sample2") == -0.375);
		REQUIRE(norm_score_diff_ratios.scores("pep3", "sample3") == -1.500);
	}
	SECTION("Test compute_ratio() works as expected")
	{
		peptide_score_data norm_score_ratio;
		norm_score_ratio.sample_names = {"sample1", "sample2", "sample3"};
		norm_score_ratio.pep_names = {"pep1", "pep2", "pep3"};

		norm_score_ratio.scores = labeled_matrix<double, std::string>(
			norm_score_ratio.pep_names.size(),
			norm_score_ratio.sample_names.size(),
			norm_score_ratio.pep_names, norm_score_ratio.sample_names
		);

		norm_score_ratio.scores("pep1", "sample1") = -8.284;
		norm_score_ratio.scores("pep1", "sample2") = 7.2849;
		norm_score_ratio.scores("pep1", "sample3") = 10.001;

		norm_score_ratio.scores("pep2", "sample1") = -0.003;
		norm_score_ratio.scores("pep2", "sample2") = 0.127482;
		norm_score_ratio.scores("pep2", "sample3") = 5.8392;

		norm_score_ratio.scores("pep3", "sample1") = 12.000;
		norm_score_ratio.scores("pep3", "sample2") = 8.89994;
		norm_score_ratio.scores("pep3", "sample3") = -4.388;

		std::unordered_map<std::string, double>pep_avgs = {
			{"pep1", 10.000}, {"pep2", 5.302}, {"pep3", -0.832}
		};

		norm_mod.compute_ratio(&norm_score_ratio, &pep_avgs);
		
		REQUIRE(norm_score_ratio.scores("pep1", "sample1") == (-8.284 / 10.000));
		REQUIRE(norm_score_ratio.scores("pep1", "sample2") == (7.2849 / 10.000));
		REQUIRE(norm_score_ratio.scores("pep1", "sample3") == (10.001 / 10.000));

		REQUIRE(norm_score_ratio.scores("pep2", "sample1") == (-0.003 / 5.302));
		REQUIRE(norm_score_ratio.scores("pep2", "sample2") == (0.127482 / 5.302));
		REQUIRE(norm_score_ratio.scores("pep2", "sample3") == (5.8392 / 5.302));

		REQUIRE(norm_score_ratio.scores("pep3", "sample1") == (12.000 / -0.832));
		REQUIRE(norm_score_ratio.scores("pep3", "sample2") == (8.89994 / -0.832));
		REQUIRE(norm_score_ratio.scores("pep3", "sample3") == (-4.388 / -0.832));
	}
}

TEST_CASE( "geometric means", "[stats]" )
{
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

TEST_CASE("Test label operations of labeled_matrix class", "[labeled_matrix]")
{
	std::vector<std::string> row_labels = {
		"Row1", "Row2", "Row3", "Row4"
	};
	std::vector<std::string> col_labels = {
		"Col1", "Col2", "Col3", "Col4"
	};
	labeled_matrix<int, std::string> labeled_mat(4, 4, row_labels, col_labels);

	SECTION("Getting row labels")
	{
		std::unordered_set<std::string> labels;
		labeled_mat.get_row_labels(labels);

		REQUIRE(labels.find(row_labels[0]) != labels.end());
		REQUIRE(labels.find(row_labels[1]) != labels.end());
		REQUIRE(labels.find(row_labels[2]) != labels.end());
		REQUIRE(labels.find(row_labels[3]) != labels.end());
	}
	SECTION("Getting column labels")
	{
		std::unordered_set<std::string> labels;
		labeled_mat.get_col_labels(labels);

		REQUIRE(labels.find(col_labels[0]) != labels.end());
		REQUIRE(labels.find(col_labels[1]) != labels.end());
		REQUIRE(labels.find(col_labels[2]) != labels.end());
		REQUIRE(labels.find(col_labels[3]) != labels.end());
	}
	SECTION("Update row labels")
	{
		std::unordered_map<std::string, std::uint32_t> labels_map;
		std::vector<std::string> labels = {
			"New Row1", "New Row2", "New Row3", "New Row4"
		};

		labeled_mat.update_row_labels(labels);
		labels_map = labeled_mat.get_row_labels();

		REQUIRE(labels_map.find(labels[0]) != labels_map.end());
		REQUIRE(labels_map.find(labels[1]) != labels_map.end());
		REQUIRE(labels_map.find(labels[2]) != labels_map.end());
		REQUIRE(labels_map.find(labels[3]) != labels_map.end());
	}
	SECTION("Update column labels")
	{
		std::unordered_map<std::string, std::uint32_t> labels_map;
		std::vector<std::string> labels = {
			"New Col1", "New Col2", "New Col3", "New Col4"
		};

		labeled_mat.update_col_labels(labels);
		labels_map = labeled_mat.get_col_labels();

		REQUIRE(labels_map.find(labels[0]) != labels_map.end());
		REQUIRE(labels_map.find(labels[1]) != labels_map.end());
		REQUIRE(labels_map.find(labels[2]) != labels_map.end());
		REQUIRE(labels_map.find(labels[3]) != labels_map.end());
	}
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

TEST_CASE("Testing swap() and filter_cols() funcitons of a labeled matrix", "[matrix], [labeled_matrix]")
{
	SECTION("Verify swapping of matrices")
	{
		matrix<double> mat = matrix<double>(3, 3);
		mat(0, 0) = 7.82;
		mat(0, 1) = 9.12;
		mat(0, 2) = 10.00;

		mat(1, 0) = -5.32;
		mat(1, 1) = 0.83;
		mat(1, 2) = -1.03;

		mat(2, 0) = 5.23;
		mat(2, 1) = 20.00;
		mat(2, 2) = 12.00;

		matrix<double> new_mat = matrix<double>(3, 4);
		new_mat(0, 0) = 0;
		new_mat(0, 1) = 0;
		new_mat(0, 2) = 0;
		new_mat(0, 3) = 0;

		new_mat(1, 0) = 0;
		new_mat(1, 1) = 0;
		new_mat(1, 2) = 0;
		new_mat(1, 3) = 0;

		new_mat(2, 0) = 0;
		new_mat(2, 1) = 0;
		new_mat(2, 2) = 0;
		new_mat(2, 3) = 0;

		swap(mat, new_mat);

		// verify mat has new_mat's values
		for (std::size_t i = 0; i < 3; i += 1)
		{
			for (std::size_t j = 0; j < 4; j += 1)
			{
				REQUIRE(mat(i, j) == 0);
			}
		}
	}
	SECTION("Verify swapping of labeled matrices")
	{
		// initialize row and column labels vectors
		std::vector<std::string> row_labels = {
			"Row1", "Row2", "Row3"
		};
		std::vector<std::string> column_labels = {
			"Column1", "Column2", "Column3"
		};
		// define a 3x3 labeled matrix with the row and column label vecs
		labeled_matrix<double, std::string>
			lab_mat = labeled_matrix<double, std::string>(
				row_labels.size(), column_labels.size(),
				row_labels, column_labels
			);
		
		// fill matrix with numbers
		lab_mat(0, 0) = -0.34;
		lab_mat(0, 1) = 8.00;
		lab_mat(0, 2) = 2.44;

		lab_mat(1, 0) = 7.01;
		lab_mat(1, 1) = 10.33;
		lab_mat(1, 2) = -3.00;
		
		lab_mat(2, 0) = 6.83;
		lab_mat(2, 1) = 8.31;
		lab_mat(2, 2) = 9.32;
		
		// initialize new row and column vectors
		std::vector<std::string> new_row_labels = {
			"New Row1", "New Row2", "New Row3"
		};
		std::vector<std::string> new_column_labels = {
			"New Column1", "New Column2", "New Column3"
		};
		// define new 3x3 labeled matrix
		labeled_matrix<double, std::string>
			new_lab_mat = labeled_matrix<double, std::string>(
				new_row_labels.size(), new_column_labels.size(),
				new_row_labels, new_column_labels
			);
		
		// fill new matrix with zeros
		new_lab_mat(0, 0) = 0;
		new_lab_mat(0, 1) = 0;
		new_lab_mat(0, 2) = 0;

		new_lab_mat(1, 0) = 0;
		new_lab_mat(1, 1) = 0;
		new_lab_mat(1, 2) = 0;

		new_lab_mat(2, 0) = 0;
		new_lab_mat(2, 1) = 0;
		new_lab_mat(2, 2) = 0;

		swap(lab_mat, new_lab_mat);

		// verify lab_mat's labels are the new labels
		std::unordered_map<std::string, std::uint32_t>
			captured_row_labels = lab_mat.get_row_labels();
		std::unordered_map<std::string, std::uint32_t>
			captured_column_labels = lab_mat.get_col_labels();
		REQUIRE(captured_row_labels.find(new_row_labels[0]) != captured_row_labels.end());
		REQUIRE(captured_row_labels.find(new_row_labels[1]) != captured_row_labels.end());
		REQUIRE(captured_row_labels.find(new_row_labels[2]) != captured_row_labels.end());

		REQUIRE(captured_column_labels.find(new_column_labels[0]) != captured_column_labels.end());
		REQUIRE(captured_column_labels.find(new_column_labels[1]) != captured_column_labels.end());
		REQUIRE(captured_column_labels.find(new_column_labels[2]) != captured_column_labels.end());

		// verify new_lab_mat's labels are the original labels
		captured_row_labels = new_lab_mat.get_row_labels();
		captured_column_labels = new_lab_mat.get_col_labels();
		REQUIRE(captured_row_labels.find(row_labels[0]) != captured_row_labels.end());
		REQUIRE(captured_row_labels.find(row_labels[1]) != captured_row_labels.end());
		REQUIRE(captured_row_labels.find(row_labels[2]) != captured_row_labels.end());

		REQUIRE(captured_column_labels.find(column_labels[0]) != captured_column_labels.end());
		REQUIRE(captured_column_labels.find(column_labels[1]) != captured_column_labels.end());
		REQUIRE(captured_column_labels.find(column_labels[2]) != captured_column_labels.end());

		// verify lab_mat has new_lab_mat's values
		for (std::size_t i = 0; i < 3; i += 1)
		{
			for (std::size_t j = 0; j < 3; j += 1)
			{
				REQUIRE(lab_mat(i, j) == 0);
			}
		}
	}
	SECTION("Verify filtering of columns")
	{
		std::vector<std::string> row_labels = {
			"Row1", "Row2", "Row3", "Row4"
		};
		std::vector<std::string> column_labels = {
			"Column1", "Column2", "Column3", "Column4"
		};
		labeled_matrix<double, std::string>
			lab_mat = labeled_matrix<double, std::string>(
				row_labels.size(), column_labels.size(),
				row_labels, column_labels
			);

		lab_mat.at(0, 0) = 0.00;
		lab_mat.at(0, 1) = 5.00;
		lab_mat.at(0, 2) = 1.00;
		lab_mat.at(0, 3) = 0.00;

		lab_mat.at(1, 0) = 1.00;
		lab_mat.at(1, 1) = 10.00;
		lab_mat.at(1, 2) = 7.00;
		lab_mat.at(1, 3) = 4.00;
		
		lab_mat.at(2, 0) = 0.00;
		lab_mat.at(2, 1) = 1.00;
		lab_mat.at(2, 2) = 2.00;
		lab_mat.at(2, 3) = 0.00;

		lab_mat.at(3, 0) = 8.00;
		lab_mat.at(3, 1) = 0.00;
		lab_mat.at(3, 2) = 1.00;
		lab_mat.at(3, 3) = 0.00;

		std::vector<std::string> filter_cols_vec = {
			"Column1", "Column4"
		};
		labeled_matrix<double, std::string> filtered_mat = lab_mat.filter_cols(filter_cols_vec);

		std::vector<std::string> captured_row_labels;
		filtered_mat.get_row_labels(captured_row_labels);
		std::vector<std::string> captured_col_labels;
		filtered_mat.get_col_labels(captured_col_labels);

		REQUIRE(captured_row_labels.size() == 4);
		REQUIRE(captured_col_labels.size() == 2);

		REQUIRE(filtered_mat.at("Row1", "Column1") == 0.00);
		REQUIRE(filtered_mat.at("Row2", "Column1") == 5.00);
		REQUIRE(filtered_mat.at("Row3", "Column1") == 1.00);
		REQUIRE(filtered_mat.at("Row4", "Column1") == 0.00);

		REQUIRE(filtered_mat.at("Row1", "Column4") == 0.00);
		REQUIRE(filtered_mat.at("Row2", "Column4") == 1.00);
		REQUIRE(filtered_mat.at("Row3", "Column4") == 10.00);
		REQUIRE(filtered_mat.at("Row4", "Column4") == 7.00);
	}
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

TEST_CASE("Relative difference", "[stats]")
{
	double a = 100;
	double b = 55;
	double expected = 0.45;
	double actual = stats::relative_difference(a, b);

	REQUIRE(actual == expected);
}

TEST_CASE("Squared difference", "[stats]")
{
	std::vector<double> double_vec = {6.0, 60.0, 10.0, 50.0};
	double subtrahend = 6.0;

	double expected = 4868.0;
	double actual = stats::squared_diff(double_vec.begin(), double_vec.end(), subtrahend);

	REQUIRE(actual == expected);
}

TEST_CASE( "Full test of Peptide Bin operations", "[peptide_bin]" )
{
    std::stringstream bins_in;
    bins_in << "pep_1\tpep_2\tpep_3\npep_4\tpep_5\tpep_6\n";
    bin_collection bin_c = peptide_bin_io::parse_bins( bins_in );

	SECTION("Test parsing and writing bins from a stream")
	{
		std::stringstream bins_out;
		peptide_bin_io::write_bins( bins_out, bin_c );

		// NOTE: test parse_bins() -> tests add_bin()
		REQUIRE( ( bin_c == peptide_bin_io::parse_bins( bins_out ) ) );
	}

	auto first_bin = bin_c.begin();
	auto second_bin = first_bin + 1;
	SECTION("Test if a bin contains a specific peptide")
	{
		REQUIRE( first_bin->contains( "pep_1" ) );
		REQUIRE( first_bin->contains( "pep_2" ) );
		REQUIRE( first_bin->contains( "pep_3" ) );
		REQUIRE( !( first_bin->contains( "pep_4" ) ) );

		REQUIRE( second_bin->contains( "pep_4" ) );
		REQUIRE( second_bin->contains( "pep_5" ) );
		REQUIRE( second_bin->contains( "pep_6" ) );
		REQUIRE( !( second_bin->contains( "pep_9" ) ) );
	}
	SECTION("Test that peptides can be added to a bin, and which bin is the smallest")
	{
		// first bin is smallest because first and second bins have same
		// number of peptides
		REQUIRE(bin_c.smallest() == *first_bin);

		first_bin->add_peptide("pep_10");
		first_bin->add_peptide("pep_11");
		first_bin->add_peptide("pep_12");

		// verify first bin is no longer smallest bin
		REQUIRE(first_bin->size() == 6);
		REQUIRE(second_bin->size() == 3);
		REQUIRE(bin_c.smallest() == *second_bin);
	}
	SECTION("Test more bins can be added, and last bin is smallest")
	{
		std::vector<peptide_bin> new_bins = {
			peptide_bin(), peptide_bin()
		};
		std::vector<std::string> other_peptides = {
			"pep24", "pep25", "pep26", "pep27"
		};

		new_bins[0].add_peptides(other_peptides.begin(), other_peptides.end());

		bin_c.add_bins(new_bins.begin(), new_bins.end());

		REQUIRE(bin_c.bins[2].size() == 4);
		REQUIRE(bin_c.bins[3].size() == 0);
		REQUIRE(bin_c.smallest() == *(bin_c.end() - 1));
	}
}

TEST_CASE("Full test of peptide scoring's individual pieces", "[peptide_scoring]")
{
	std::vector<double> expected_scores = {
		0.00, 10.00, 8.00, 15.00, 3.00, 1.00, 3.00, 
		0.00, 10.00, 5.00, 8.00, 2.00, 0.00, 2.00, 
		21.00, 27.00, 9.00, 24.00, 14.00, 0.00, 14.00, 
		2.00, 6.00, 7.00, 7.00, 1.00, 0.00, 5.00, 
		2.00, 10.00, 31.00, 22.00, 1.00, 1.00, 14.00, 
		12.00, 28.00, 24.00, 19.00, 13.00, 0.00, 23.00, 
		0.00, 14.00, 11.00, 18.00, 6.00, 1.00, 5.00, 
		2.00, 17.00, 19.00, 23.00, 4.00, 2.00, 6.00, 
		2.00, 13.00, 11.00, 25.00, 0.00, 1.00, 5.00, 
		1.00, 24.00, 15.00, 30.00, 9.00, 1.00, 10.00
	};
	peptide_score_data score_data;
	peptide_scoring::parse_peptide_scores(
		score_data,
		"../test/input_data/test_enrich_raw_scores.tsv"
	);

	SECTION("Verify file name of score data points to test/test_enrich_raw_scores.tsv")
	{
		REQUIRE(score_data.file_name.compare("test_enrich_raw_scores.tsv") == 0);
	}
	SECTION("Verify peptide names of score data match those found in raw scores file")
	{
		std::vector<std::string> expected_pep_names = {
			"PV1_000000", "PV1_000001", "PV1_000002", "PV1_000003", "PV1_000004", "PV1_000005", "PV1_000006"
		};

		REQUIRE(expected_pep_names.size() == score_data.pep_names.size());
		for (std::size_t i = 0; i < expected_pep_names.size(); i += 1)
		{
			REQUIRE(expected_pep_names[i].compare(score_data.pep_names[i]) == 0);
		}
	}
	SECTION("Verify sample names of score data match those found in raw scores file")
	{
		std::vector<std::string> expected_sample_names = {
			"HCV-I_PV1_Run3_B", "HCV-H_PV1_Run3_B", "HCV-J_PV1_Run3_A",
			"JAS_PV1_Run10_A", "HCV-H_PV1_Run3_A", "HCV-I_PV1_Run3_A",
			"HCV-J_PV1_Run3_B", "SB_PV1_Run3_A", "JAS_PV1_Run10_B",
			"SB_PV1_Run3_B"
		};

		REQUIRE(expected_sample_names.size() == score_data.sample_names.size());
		for (std::size_t i = 0; i < expected_sample_names.size(); i += 1)
		{
			REQUIRE(expected_sample_names[i].compare(score_data.sample_names[i]) == 0);
		}
	}
	SECTION("Verify scores match those found in raw scores file")
	{
		for (std::size_t i = 0; i < score_data.sample_names.size(); i += 1)
		{
			for (std::size_t j = 0; j < score_data.pep_names.size(); j += 1)
			{
				REQUIRE(score_data.scores(i, j) == expected_scores[i * score_data.pep_names.size() + j]);
			}
		}
	}
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
			// verify getting probe(s) by rank
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
			// verify getting probe(s) by rank
			REQUIRE(*int_probe_rank.get_probes_with_rank(1.0) == *expected_values.find(1.0));
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
			// verify getting probe(s) by rank
			REQUIRE(*int_probe_rank.get_probes_with_rank(1.05) == *expected_values.find(1.05));
			REQUIRE(*int_probe_rank.get_probes_with_rank(156.10) == *expected_values.find(156.10));
			REQUIRE(*int_probe_rank.get_probes_with_rank(1.4522) == *expected_values.find(1.45));
        }
}

TEST_CASE("Test of get_raw_sums() and get_raw_scores()", "[module_enrich]")
{
	module_enrich enrich_mod;
	SECTION("Verifying collection of raw scores")
	{
        std::vector<std::vector<double>> raw_scores_vec;
        peptide_score_data_sample_major raw_scores_dat;
        peptide_score_data_sample_major* raw_scores_dat_ptr = nullptr;
        peptide_scoring::parse_peptide_scores(
            raw_scores_dat, "../test/input_data/test_enrich_raw_scores.tsv"
        );
        raw_scores_dat_ptr = &raw_scores_dat;
        raw_scores_vec.resize(raw_scores_dat.scores.ncols());

        std::vector<std::string> sample_names = {
            "HCV-I_PV1_Run3_B", "JAS_PV1_Run10_A", "HCV-J_PV1_Run3_B"
        };

        enrich_mod.get_raw_scores(&raw_scores_vec, raw_scores_dat_ptr, sample_names);

        REQUIRE(raw_scores_vec[0][0] == 0.00);
        REQUIRE(raw_scores_vec[0][1] == 2.00);
        REQUIRE(raw_scores_vec[0][2] == 0.00);

        REQUIRE(raw_scores_vec[1][0] == 10.00);
        REQUIRE(raw_scores_vec[1][1] == 6.00);
        REQUIRE(raw_scores_vec[1][2] == 14.00);

        REQUIRE(raw_scores_vec[2][0] == 8.00);
        REQUIRE(raw_scores_vec[2][1] == 7.00);
        REQUIRE(raw_scores_vec[2][2] == 11.00);

        REQUIRE(raw_scores_vec[3][0] == 15.00);
        REQUIRE(raw_scores_vec[3][1] == 7.00);
        REQUIRE(raw_scores_vec[3][2] == 18.00);

        REQUIRE(raw_scores_vec[4][0] == 3.00);
        REQUIRE(raw_scores_vec[4][1] == 1.00);
        REQUIRE(raw_scores_vec[4][2] == 6.00);

        REQUIRE(raw_scores_vec[5][0] == 1.00);
        REQUIRE(raw_scores_vec[5][1] == 0.00);
        REQUIRE(raw_scores_vec[5][2] == 1.00);

        REQUIRE(raw_scores_vec[6][0] == 3.00);
        REQUIRE(raw_scores_vec[6][1] == 5.00);
        REQUIRE(raw_scores_vec[6][2] == 5.00);
	}
	SECTION("Verifying collection of raw sums")
	{
        std::vector<double> col_sums;
        std::vector<std::vector<double>> raw_scores_vec = {
            {0, 0, 21, 2, 2, 12, 0, 2, 2, 1},
            {10, 10, 27, 6, 10, 28, 14, 17, 13, 24},
            {8, 5, 9, 7, 31, 24, 11, 19, 11, 15},
            {15, 8, 24, 7, 22, 19, 18, 23, 25, 30},
            {3, 2, 14, 1, 1, 13, 6, 4, 0, 9},
            {1, 0, 0, 0, 1, 0, 1, 2, 1, 1},
            {3, 2, 14, 5, 14, 23, 5, 6, 5, 10}
        };

        col_sums = enrich_mod.get_raw_sums(raw_scores_vec);

        REQUIRE(col_sums[0] == 40.00);
        REQUIRE(col_sums[1] == 27.00);
        REQUIRE(col_sums[2] == 109.00);
        REQUIRE(col_sums[3] == 28.00);
        REQUIRE(col_sums[4] == 81.00);
        REQUIRE(col_sums[5] == 119.00);
        REQUIRE(col_sums[6] == 55.00);
        REQUIRE(col_sums[7] == 73.00);
        REQUIRE(col_sums[8] == 57.00);
        REQUIRE(col_sums[9] == 90.00);
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

/* remove
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
remove */

TEST_CASE("Full test of subjoin's individual methods", "[module_subjoin]")
{
	module_subjoin sub_mod;
	SECTION("Test of namelist parsing")
	{
		// initialize in-file stream with path to namelist
		std::ifstream namelist("../test/input_data/test_subjoin_namelist.txt", std::ios_base::in);
		// define destination vector of strings
		std::vector<std::string> names;
		// name replacement map
		module_subjoin::name_replacement_list replacement_map;

		replacement_map = sub_mod.parse_namelist(names, namelist);

		REQUIRE(names[0].compare("GY_PV_Run3X") == 0);
		REQUIRE(names[1].compare("GY_PV_Run4X") == 0);
		REQUIRE(names[2].compare("VKi_6vk_A") == 0);
		REQUIRE(names[3].compare("VKi_6vk_B") == 0);
		REQUIRE(names[4].compare("URB-PV1_SUB_6_A") == 0);
		REQUIRE(names[5].compare("URB-PV1_SUB_6_B") == 0);
		REQUIRE(names[6].compare("URB-PV1_SUB_6_C") == 0);

		for (
			std::size_t name_idx = 0;
			name_idx < replacement_map.size();
			name_idx += 1
		)
		{
			REQUIRE(replacement_map.find(names[name_idx]) != replacement_map.end());
		}
	}

	// initialize first peptide score sample major with data
	peptide_score_data first_data_set;
	first_data_set.file_name = "first_test_data_set";
	first_data_set.pep_names = {"PEP_006803", "PEP_000039", "PEP_003001"};
	first_data_set.sample_names = {"RX-PV1_06", "RX-PV1_02", "RX-PV1_03"};

	first_data_set.scores = labeled_matrix<double, std::string>(
		first_data_set.pep_names.size(), first_data_set.sample_names.size(),
		first_data_set.pep_names, first_data_set.sample_names
	);

	first_data_set.scores(0, 0) = 0.00;
	first_data_set.scores(0, 1) = 5.00;
	first_data_set.scores(0, 2) = 25.00;

	first_data_set.scores(1, 0) = 10.00;
	first_data_set.scores(1, 1) = 30.00;
	first_data_set.scores(1, 2) = 112.00;

	first_data_set.scores(2, 0) = 18.00;
	first_data_set.scores(2, 1) = 12.00;
	first_data_set.scores(2, 2) = 4.00;

	// initialize second peptide score sample major with data
	peptide_score_data second_data_set;
	second_data_set.file_name = "second_test_data_set";
	second_data_set.pep_names = {"PEP_006803", "PEP_000039", "PEP_003001"};
	second_data_set.sample_names = {"RX-PV1_06", "RX-PV1_12", "RX-PV1_04"};
	
	second_data_set.scores = labeled_matrix<double, std::string>(
		second_data_set.pep_names.size(),
		second_data_set.sample_names.size(),
		second_data_set.pep_names, second_data_set.sample_names
	);

	second_data_set.scores(0, 0) = 14.00;
	second_data_set.scores(0, 1) = 10.00;
	second_data_set.scores(0, 2) = 10.00;

	second_data_set.scores(1, 0) = 12.00;
	second_data_set.scores(1, 1) = 9.00;
	second_data_set.scores(1, 2) = 20.00;

	second_data_set.scores(2, 0) = 0.00;
	second_data_set.scores(2, 1) = 57.00;
	second_data_set.scores(2, 2) = 35.00;

	SECTION("Test of joining with INGORE resolution strategy")
	{
		// initialize IGNORE score strategy
		options_subjoin s_opts;
		s_opts.duplicate_resolution_strategy = evaluation_strategy::duplicate_resolution_strategy::IGNORE;
	
		// define destination labeled matrix associated doubles with strings
		labeled_matrix<double, std::string> joined_matrix;

		joined_matrix = sub_mod.join_with_resolve_strategy(
			first_data_set, second_data_set,
			s_opts.duplicate_resolution_strategy
		);

		std::vector<std::string> row_labels;
		joined_matrix.get_row_labels(row_labels);
		std::vector<std::string> col_labels;
		joined_matrix.get_col_labels(col_labels);

		REQUIRE(row_labels.size() == 3);
		REQUIRE(col_labels.size() == 5);

		REQUIRE(joined_matrix("PEP_003001", "RX-PV1_03") == 4.00);
		REQUIRE(joined_matrix("PEP_000039", "RX-PV1_03") == 112.00);
		REQUIRE(joined_matrix("PEP_006803", "RX-PV1_03") == 25.00);

		REQUIRE(joined_matrix("PEP_003001", "RX-PV1_04") == 35.00);
		REQUIRE(joined_matrix("PEP_000039", "RX-PV1_04") == 20.00);
		REQUIRE(joined_matrix("PEP_006803", "RX-PV1_04") == 10.00);

		REQUIRE(joined_matrix("PEP_003001", "RX-PV1_12") == 57.00);
		REQUIRE(joined_matrix("PEP_000039", "RX-PV1_12") == 9.00);
		REQUIRE(joined_matrix("PEP_006803", "RX-PV1_12") == 10.00);
		
		REQUIRE(joined_matrix("PEP_003001", "RX-PV1_06") == 0.00);
		REQUIRE(joined_matrix("PEP_000039", "RX-PV1_06") == 12.00);
		REQUIRE(joined_matrix("PEP_006803", "RX-PV1_06") == 14.00);

		REQUIRE(joined_matrix("PEP_003001", "RX-PV1_02") == 12.00);
		REQUIRE(joined_matrix("PEP_000039", "RX-PV1_02") == 30.00);
		REQUIRE(joined_matrix("PEP_006803", "RX-PV1_02") == 5.00);
	}
	SECTION("Test of joining with COMBINE resolution strategy")
	{
		// initialize COMBINE score strategy
		options_subjoin s_opts;
		s_opts.duplicate_resolution_strategy = evaluation_strategy::duplicate_resolution_strategy::COMBINE;

		// define destination labeled matrix associated doubles with strings
		labeled_matrix<double, std::string> joined_matrix;

		joined_matrix = sub_mod.join_with_resolve_strategy(
			first_data_set, second_data_set,
			s_opts.duplicate_resolution_strategy
		);

		std::vector<std::string> row_labels;
		joined_matrix.get_row_labels(row_labels);
		std::vector<std::string> col_labels;
		joined_matrix.get_col_labels(col_labels);

		REQUIRE(row_labels.size() == 3);
		REQUIRE(col_labels.size() == 5);

		REQUIRE(joined_matrix("PEP_003001", "RX-PV1_03") == 4.00);
		REQUIRE(joined_matrix("PEP_000039", "RX-PV1_03") == 112.00);
		REQUIRE(joined_matrix("PEP_006803", "RX-PV1_03") == 25.00);

		REQUIRE(joined_matrix("PEP_003001", "RX-PV1_04") == 35.00);
		REQUIRE(joined_matrix("PEP_000039", "RX-PV1_04") == 20.00);
		REQUIRE(joined_matrix("PEP_006803", "RX-PV1_04") == 10.00);

		REQUIRE(joined_matrix("PEP_003001", "RX-PV1_12") == 57.00);
		REQUIRE(joined_matrix("PEP_000039", "RX-PV1_12") == 9.00);
		REQUIRE(joined_matrix("PEP_006803", "RX-PV1_12") == 10.00);

		REQUIRE(joined_matrix("PEP_003001", "RX-PV1_06") == 18.00);
		REQUIRE(joined_matrix("PEP_000039", "RX-PV1_06") == 22.00);
		REQUIRE(joined_matrix("PEP_006803", "RX-PV1_06") == 14.00);

		REQUIRE(joined_matrix("PEP_003001", "RX-PV1_02") == 12.00);
		REQUIRE(joined_matrix("PEP_000039", "RX-PV1_02") == 30.00);
		REQUIRE(joined_matrix("PEP_006803", "RX-PV1_02") == 5.00);
	}
	SECTION("Test of joining with INCLUDE resolution strategy")
	{
		// initialize INCLUDE score strategy
		options_subjoin s_opts;
		s_opts.duplicate_resolution_strategy = evaluation_strategy::duplicate_resolution_strategy::INCLUDE;

		// define destination labeled matrix associated doubles with strings
		labeled_matrix<double, std::string> joined_matrix;

		joined_matrix = sub_mod.join_with_resolve_strategy(
			first_data_set, second_data_set,
			s_opts.duplicate_resolution_strategy
		);

		std::vector<std::string> row_labels;
		joined_matrix.get_row_labels(row_labels);
		std::vector<std::string> col_labels;
		joined_matrix.get_col_labels(col_labels);

		REQUIRE(row_labels.size() == 6);
		REQUIRE(col_labels.size() == 5);

		REQUIRE(joined_matrix("PEP_003001_first_test_data_set", "RX-PV1_03") == 4.00);
		REQUIRE(joined_matrix("PEP_006803_second_test_data_set", "RX-PV1_03") == 0.00);
		REQUIRE(joined_matrix("PEP_003001_second_test_data_set", "RX-PV1_03") == 0.00);
		REQUIRE(joined_matrix("PEP_000039_second_test_data_set", "RX-PV1_03") == 0.00);
		REQUIRE(joined_matrix("PEP_000039_first_test_data_set", "RX-PV1_03") == 112.00);
		REQUIRE(joined_matrix("PEP_006803_first_test_data_set", "RX-PV1_03") == 25.00);

		REQUIRE(joined_matrix("PEP_003001_first_test_data_set", "RX-PV1_04") == 0.00);
		REQUIRE(joined_matrix("PEP_006803_second_test_data_set", "RX-PV1_04") == 10.00);
		REQUIRE(joined_matrix("PEP_003001_second_test_data_set", "RX-PV1_04") == 35.00);
		REQUIRE(joined_matrix("PEP_000039_second_test_data_set", "RX-PV1_04") == 20.00);
		REQUIRE(joined_matrix("PEP_000039_first_test_data_set", "RX-PV1_04") == 0.00);
		REQUIRE(joined_matrix("PEP_006803_first_test_data_set", "RX-PV1_04") == 0.00);

		REQUIRE(joined_matrix("PEP_003001_first_test_data_set", "RX-PV1_12") == 0.00);
		REQUIRE(joined_matrix("PEP_006803_second_test_data_set", "RX-PV1_12") == 10.00);
		REQUIRE(joined_matrix("PEP_003001_second_test_data_set", "RX-PV1_12") == 57.00);
		REQUIRE(joined_matrix("PEP_000039_second_test_data_set", "RX-PV1_12") == 9.00);
		REQUIRE(joined_matrix("PEP_000039_first_test_data_set", "RX-PV1_12") == 0.00);
		REQUIRE(joined_matrix("PEP_006803_first_test_data_set", "RX-PV1_12") == 0.00);

		REQUIRE(joined_matrix("PEP_003001_first_test_data_set", "RX-PV1_06") == 18.00);
		REQUIRE(joined_matrix("PEP_006803_second_test_data_set", "RX-PV1_06") == 14.00);
		REQUIRE(joined_matrix("PEP_003001_second_test_data_set", "RX-PV1_06") == 0.00);
		REQUIRE(joined_matrix("PEP_000039_second_test_data_set", "RX-PV1_06") == 12.00);
		REQUIRE(joined_matrix("PEP_000039_first_test_data_set", "RX-PV1_06") == 10.00);
		REQUIRE(joined_matrix("PEP_006803_first_test_data_set", "RX-PV1_06") == 0.00);

		REQUIRE(joined_matrix("PEP_003001_first_test_data_set", "RX-PV1_02") == 12.00);
		REQUIRE(joined_matrix("PEP_006803_second_test_data_set", "RX-PV1_02") == 0.00);
		REQUIRE(joined_matrix("PEP_003001_second_test_data_set", "RX-PV1_02") == 0.00);
		REQUIRE(joined_matrix("PEP_000039_second_test_data_set", "RX-PV1_02") == 0.00);
		REQUIRE(joined_matrix("PEP_000039_first_test_data_set", "RX-PV1_02") == 30.00);
		REQUIRE(joined_matrix("PEP_006803_first_test_data_set", "RX-PV1_02") == 5.00);
	}
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

TEST_CASE( "Run Subjoin exclude option", "[module_subjoin]" )
{
    module_subjoin mod;
    options_subjoin opts;
    opts.exclude_names = true;
    opts.use_sample_names = true;
    opts.out_matrix_fname = "../test/test_subjoin_exclude_output.tsv";
    opts.input_matrix_name_pairs.emplace_back( std::make_pair( "../test/input_data/test_zscore_score_matrix.tsv", 
    															"../test/input_data/test_subjoin_exclude_namelist.txt" ) );
    mod.run( &opts );

    std::string expected = "../test/expected/test_expected_subjoin_exclude_output.tsv";
    std::string actual = "../test/test_subjoin_exclude_output.tsv";
    std::ifstream ifexpected( expected, std::ios_base::in );
    std::ifstream ifactual( actual, std::ios_base::in );
    std::string expected_line;
    std::string actual_line;
    std::unordered_set<std::string> expected_set;
    std::unordered_set<std::string> actual_set;

    // add each to the set
	while( std::getline(ifexpected, expected_line) && std::getline(ifactual, actual_line) )
        {	
            expected_set.insert( expected_line );
            actual_set.insert( actual_line );
        }
	ifexpected.close();
	ifactual.close();

	REQUIRE( expected_set == actual_set );
}

TEST_CASE("Verify metadata map construction operation", "[metadata_map]")
{
    std::string metadata_map_fname = "../test/input_data/full_design_clean_min30_taxtweak_100perc_jingmens_2019-09-12_segment.metadata";
}

TEST_CASE("Full test of link module's individual methods", "[module_link]")
{
	module_link link_mod;

	SECTION("Testing prot map creation")
	{
        std::vector<std::vector<scored_entity<std::string, double>>> expected_entities = {
            { // GTGA
                scored_entity<std::string, double>("2232", 1.00)
            },
            { // AGGA
                scored_entity<std::string, double>("4857", 1.00)
            },
            { // GTTA
                scored_entity<std::string, double>("1432", 1.00)
            },
            { // AAAG
                scored_entity<std::string, double>("8453", 1.00),
                scored_entity<std::string, double>("1432", 1.00),
                scored_entity<std::string, double>("4857", 1.00),
                scored_entity<std::string, double>("2232", 1.00)
            },
            { // AGTT
                scored_entity<std::string, double>("1432", 1.00),
                scored_entity<std::string, double>("2232", 1.00)
            },
            { // TTAG
                scored_entity<std::string, double>("1432", 1.00),
                scored_entity<std::string, double>("4857", 1.00)
            },
            { // AAGT
                scored_entity<std::string, double>("1432", 1.00),
                scored_entity<std::string, double>("2232", 1.00)
            },
            { // TAGG
                scored_entity<std::string, double>("4857", 1.00)
            },
            { // GGAA
                scored_entity<std::string, double>("4857", 1.00)
            },
            { // GAAG
                scored_entity<std::string, double>("8453", 1.00)
            },
            { // TGAA
                scored_entity<std::string, double>("8453", 1.00),
                scored_entity<std::string, double>("2232", 1.00)
            },
            { // GAAA
                scored_entity<std::string, double>("8453", 1.00),
                scored_entity<std::string, double>("1432", 1.00),
                scored_entity<std::string, double>("4857", 1.00),
                scored_entity<std::string, double>("2232", 1.00)
            },
            { // AAGA
                scored_entity<std::string, double>("8453", 1.00)
            },
            { // AGAA
                scored_entity<std::string, double>("8453", 1.00)
            }
        };
        std::unordered_map<
            std::string,
            std::unordered_set<scored_entity<std::string, double>>
        > dest_kmer_map;

		std::vector<sequence> seq_vec = {
			sequence(
				">ID=AHGTV_HH1 AC=ATFGDDVAHS_6TT OXX=4672,5934,8891,8453",
				"TGAAGAAAG"
			),
			sequence(
				">ID=AVFTGS_RR3 AC=UJDNNDFMT_PI6 OXX=5864,1342,0795,1432",
				"GAAAGTTAG"
			),
			sequence(
				">ID=GTHYVS_PPE AC=GWVBHW5FH_WT6 OXX=7284,8492,1189,4857",
				"TTAGGAAAG"
			),
			sequence(
				">ID=XCCSBV_THS AC=QQWVTH6HY_EHR OXX=5362,0958,4782,2232",
				"GTGAAAGTT"
			)
		};

		link_mod.create_prot_map<std::size_t>(
			dest_kmer_map, seq_vec, 4, 3
		);

        REQUIRE(dest_kmer_map.size() == 14);

        auto found = dest_kmer_map.find("GTGA");
        REQUIRE(found != dest_kmer_map.end());
        auto entity_set = found->second;
        REQUIRE(entity_set.find(expected_entities[0][0]) != entity_set.end());

        found = dest_kmer_map.find("AGGA");
        REQUIRE(found != dest_kmer_map.end());
        entity_set = found->second;
        REQUIRE(entity_set.find(expected_entities[1][0]) != entity_set.end());

        found = dest_kmer_map.find("GTTA");
        REQUIRE(found != dest_kmer_map.end());
        entity_set = found->second;
        REQUIRE(entity_set.find(expected_entities[2][0]) != entity_set.end());

        found = dest_kmer_map.find("AAAG");
        REQUIRE(found != dest_kmer_map.end());
        entity_set = found->second;
        REQUIRE(entity_set.find(expected_entities[3][0]) != entity_set.end());
        REQUIRE(entity_set.find(expected_entities[3][1]) != entity_set.end());
        REQUIRE(entity_set.find(expected_entities[3][2]) != entity_set.end());
        REQUIRE(entity_set.find(expected_entities[3][3]) != entity_set.end());

        found = dest_kmer_map.find("AGTT");
        REQUIRE(found != dest_kmer_map.end());
        entity_set = found->second;
        REQUIRE(entity_set.find(expected_entities[4][0]) != entity_set.end());
        REQUIRE(entity_set.find(expected_entities[4][1]) != entity_set.end());

        found = dest_kmer_map.find("TTAG");
        REQUIRE(found != dest_kmer_map.end());
        entity_set = found->second;
        REQUIRE(entity_set.find(expected_entities[5][0]) != entity_set.end());
        REQUIRE(entity_set.find(expected_entities[5][1]) != entity_set.end());

        found = dest_kmer_map.find("AAGT");
        REQUIRE(found != dest_kmer_map.end());
        entity_set = found->second;
        REQUIRE(entity_set.find(expected_entities[6][0]) != entity_set.end());
        REQUIRE(entity_set.find(expected_entities[6][1]) != entity_set.end());

        found = dest_kmer_map.find("TAGG");
        REQUIRE(found != dest_kmer_map.end());
        entity_set = found->second;
        REQUIRE(entity_set.find(expected_entities[7][0]) != entity_set.end());

        found = dest_kmer_map.find("GGAA");
        REQUIRE(found != dest_kmer_map.end());
        entity_set = found->second;
        REQUIRE(entity_set.find(expected_entities[8][0]) != entity_set.end());

        found = dest_kmer_map.find("GAAG");
        REQUIRE(found != dest_kmer_map.end());
        entity_set = found->second;
        REQUIRE(entity_set.find(expected_entities[9][0]) != entity_set.end());

        found = dest_kmer_map.find("TGAA");
        REQUIRE(found != dest_kmer_map.end());
        entity_set = found->second;
        REQUIRE(entity_set.find(expected_entities[10][0]) != entity_set.end());
        REQUIRE(entity_set.find(expected_entities[10][1]) != entity_set.end());

        found = dest_kmer_map.find("GAAA");
        REQUIRE(found != dest_kmer_map.end());
        entity_set = found->second;
        REQUIRE(entity_set.find(expected_entities[11][0]) != entity_set.end());
        REQUIRE(entity_set.find(expected_entities[11][1]) != entity_set.end());
        REQUIRE(entity_set.find(expected_entities[11][2]) != entity_set.end());
        REQUIRE(entity_set.find(expected_entities[11][3]) != entity_set.end());

        found = dest_kmer_map.find("AAGA");
        REQUIRE(found != dest_kmer_map.end());
        entity_set = found->second;
        REQUIRE(entity_set.find(expected_entities[12][0]) != entity_set.end());

        found = dest_kmer_map.find("AGAA");
        REQUIRE(found != dest_kmer_map.end());
        entity_set = found->second;
        REQUIRE(entity_set.find(expected_entities[13][0]) != entity_set.end());
	}
	SECTION("Testing pep map creation")
	{
        std::unordered_map<
            std::string,
            std::unordered_set<scored_entity<std::string, double>>
        > kmer_map = {
            {
                "AAAG",
                {
                    scored_entity<std::string, double>("8453", 1.00), scored_entity<std::string, double>("1432", 1.00),
                    scored_entity<std::string, double>("4857", 1.00), scored_entity<std::string, double>("2232", 1.00)
                }
            },
            {
                "TGAA",
                {
                    scored_entity<std::string, double>("8453", 1.00), scored_entity<std::string, double>("4857", 1.00),
                    scored_entity<std::string, double>("2232", 1.00)
                }
            }
        };

		std::vector<sequence> seq_vec = {
			sequence(
				">ID=AHGTV_HH1 AC=ATFGDDVAHS_6TT OXX=4672,5934,8891,8453",
				"ATTAACCTGGTCGTGAAAG"
			),
			sequence(
				">ID=AVFTGS_RR3 AC=UJDNNDFMT_PI6 OXX=5864,1342,0795,1432",
				"GTGAAAGGGTCCCTGGAAAG"
			),
			sequence(
				">ID=GTHYVS_PPE AC=GWVBHW5FH_WT6 OXX=7284,8492,1189,4857",
				"TTAGATTGTGTCCCGAAAGT"
			),
			sequence(
				">ID=XCCSBV_THS AC=QQWVTH6HY_EHR OXX=5362,0958,4782,2232",
				"TGGTCGTGAAAGTTAGATTA"
			)
		};

        std::vector<std::vector<scored_entity<std::string, double>>> expected_entities = {
            {
                scored_entity<std::string, double>("2232", 2.00),
                scored_entity<std::string, double>("4857", 2.00),
                scored_entity<std::string, double>("8453", 2.00),
                scored_entity<std::string, double>("1432", 1.00)
            },
            {
                scored_entity<std::string, double>("2232", 3.00),
                scored_entity<std::string, double>("4857", 3.00),
                scored_entity<std::string, double>("8453", 3.00),
                scored_entity<std::string, double>("1432", 2.00)
            },
            {
                scored_entity<std::string, double>("2232", 1.00),
                scored_entity<std::string, double>("4857", 1.00),
                scored_entity<std::string, double>("8453", 1.00),
                scored_entity<std::string, double>("1432", 1.00)
            },
            {
                scored_entity<std::string, double>("2232", 2.00),
                scored_entity<std::string, double>("4857", 2.00),
                scored_entity<std::string, double>("8453", 2.00),
                scored_entity<std::string, double>("1432", 1.00)
            }
        };

		std::vector<
			std::tuple<
				std::string,
                std::unordered_set<
				    scored_entity<std::string, double>
				>
			>
		> dest_pep_vec;
		link_mod.create_pep_map(kmer_map, dest_pep_vec, seq_vec, 4);

        REQUIRE(dest_pep_vec.size() == 4);

        auto pep = dest_pep_vec.begin();
        auto entity_set = std::get<1>(*pep);
        REQUIRE(std::get<0>(*pep).compare(">ID=AHGTV_HH1 AC=ATFGDDVAHS_6TT OXX=4672,5934,8891,8453") == 0);
        REQUIRE(entity_set.find(expected_entities[0][0]) != entity_set.end());
        REQUIRE(entity_set.find(expected_entities[0][1]) != entity_set.end());
        REQUIRE(entity_set.find(expected_entities[0][2]) != entity_set.end());
        REQUIRE(entity_set.find(expected_entities[0][3]) != entity_set.end());

        pep++;
        entity_set = std::get<1>(*pep);
        REQUIRE(std::get<0>(*pep).compare(">ID=AVFTGS_RR3 AC=UJDNNDFMT_PI6 OXX=5864,1342,0795,1432") == 0);
        REQUIRE(entity_set.find(expected_entities[1][0]) != entity_set.end());
        REQUIRE(entity_set.find(expected_entities[1][1]) != entity_set.end());
        REQUIRE(entity_set.find(expected_entities[1][2]) != entity_set.end());
        REQUIRE(entity_set.find(expected_entities[1][3]) != entity_set.end());

        pep++;
        entity_set = std::get<1>(*pep);
        REQUIRE(std::get<0>(*pep).compare(">ID=GTHYVS_PPE AC=GWVBHW5FH_WT6 OXX=7284,8492,1189,4857") == 0);
        REQUIRE(entity_set.find(expected_entities[2][0]) != entity_set.end());
        REQUIRE(entity_set.find(expected_entities[2][1]) != entity_set.end());
        REQUIRE(entity_set.find(expected_entities[2][2]) != entity_set.end());
        REQUIRE(entity_set.find(expected_entities[2][3]) != entity_set.end());

        pep++;
        entity_set = std::get<1>(*pep);
        REQUIRE(std::get<0>(*pep).compare(">ID=XCCSBV_THS AC=QQWVTH6HY_EHR OXX=5362,0958,4782,2232") == 0);
        REQUIRE(entity_set.find(expected_entities[3][0]) != entity_set.end());
        REQUIRE(entity_set.find(expected_entities[3][1]) != entity_set.end());
        REQUIRE(entity_set.find(expected_entities[3][2]) != entity_set.end());
        REQUIRE(entity_set.find(expected_entities[3][3]) != entity_set.end());
	}
	SECTION("Testing pep map creation with kmer penalty")
	{
        std::unordered_map<
            std::string,
            std::unordered_set<scored_entity<std::string, double>>
        > kmer_map = {
            {
                "AAAG",
                {
                    scored_entity<std::string, double>("8453", 1.00),
                    scored_entity<std::string, double>("1432", 1.00),
                    scored_entity<std::string, double>("4857", 1.00),
                    scored_entity<std::string, double>("2232", 1.00)
                }
            },
            {
                "TGAA",
                {
                    scored_entity<std::string, double>("8453", 1.00),
                    scored_entity<std::string, double>("4857", 1.00),
                    scored_entity<std::string, double>("2232", 1.00)
                }
            }
        };

        std::vector<sequence> seq_vec = {
            sequence(
                ">ID=XFYYS_THS AC=YYWHC6_TG OXX=7823,1048,2938",
                "ATTGACTGAA"
            ),
            sequence(
                ">ID=VVBT03_23 AC=BAUCH77_T OXX=1242,9991,8991",
                "GAAAGTGAAC"
            ),
            sequence(
                ">ID=BNYSS_06S AC=HISTAA_WQ OXX=8923,1213,0019",
                "GATGAACGAT"
            ),
            sequence(
                ">ID=HHNIS_01S AC=UUNOUU_12 OXX=9901,1248,2312",
                "GTACAAAAAG"
            )
        };

        std::vector<std::vector<scored_entity<std::string, double>>>
            expected_entities = {
                {
                    scored_entity<std::string, double>("8453", 0.33),
                    scored_entity<std::string, double>("4857", 0.33),
                    scored_entity<std::string, double>("2232", 0.33)
                },
                {
                    scored_entity<std::string, double>("8453", 0.83),
                    scored_entity<std::string, double>("1432", 0.50),
                    scored_entity<std::string, double>("4857", 0.83),
                    scored_entity<std::string, double>("2232", 0.83)
                },
                {
                    scored_entity<std::string, double>("8453", 0.33),
                    scored_entity<std::string, double>("4857", 0.33),
                    scored_entity<std::string, double>("2232", 0.33)
                },
                {
                    scored_entity<std::string, double>("8453", 0.50),
                    scored_entity<std::string, double>("1432", 0.50),
                    scored_entity<std::string, double>("4857", 0.50),
                    scored_entity<std::string, double>("2232", 0.50)
                }
            };

        std::vector<
            std::tuple<
                std::string,
                std::unordered_set<
                    scored_entity<std::string, double>
                >
            >
        > dest_pep_vec;
        link_mod.create_pep_map_with_kmer_penalty(kmer_map, dest_pep_vec, seq_vec, 4);

        REQUIRE(dest_pep_vec.size() == 4);

        auto pep = dest_pep_vec.begin();
        auto entity_set = std::get<1>(*pep);
        REQUIRE(std::get<0>(*pep).compare(">ID=XFYYS_THS AC=YYWHC6_TG OXX=7823,1048,2938") == 0);
        REQUIRE(entity_set.find(expected_entities[0][0]) != entity_set.end());
        REQUIRE(entity_set.find(expected_entities[0][1]) != entity_set.end());
        REQUIRE(entity_set.find(expected_entities[0][2]) != entity_set.end());

        pep++;
        entity_set = std::get<1>(*pep);
        REQUIRE(std::get<0>(*pep).compare(">ID=VVBT03_23 AC=BAUCH77_T OXX=1242,9991,8991") == 0);
        REQUIRE(entity_set.find(expected_entities[1][0]) != entity_set.end());
        REQUIRE(entity_set.find(expected_entities[1][1]) != entity_set.end());
        REQUIRE(entity_set.find(expected_entities[1][2]) != entity_set.end());
        REQUIRE(entity_set.find(expected_entities[1][3]) != entity_set.end());

        pep++;
        entity_set = std::get<1>(*pep);
        REQUIRE(std::get<0>(*pep).compare(">ID=BNYSS_06S AC=HISTAA_WQ OXX=8923,1213,0019") == 0);
        REQUIRE(entity_set.find(expected_entities[2][0]) != entity_set.end());
        REQUIRE(entity_set.find(expected_entities[2][1]) != entity_set.end());
        REQUIRE(entity_set.find(expected_entities[2][2]) != entity_set.end());

        pep++;
        entity_set = std::get<1>(*pep);
        REQUIRE(std::get<0>(*pep).compare(">ID=HHNIS_01S AC=UUNOUU_12 OXX=9901,1248,2312") == 0);
        REQUIRE(entity_set.find(expected_entities[3][0]) != entity_set.end());
        REQUIRE(entity_set.find(expected_entities[3][1]) != entity_set.end());
        REQUIRE(entity_set.find(expected_entities[3][2]) != entity_set.end());
        REQUIRE(entity_set.find(expected_entities[3][3]) != entity_set.end());
	}
	SECTION("Testing verification of ID type with ID index")
	{
        // subverts ambiguity when calling verify_id_type() with 0
        auto verify_id_wrapper = [&link_mod](std::string name, std::size_t retr)
        {
            return link_mod.verify_id_type(name, retr);
        };

		sequence seq(
			">ID=HNAC_EYY AC=KANMEA_R45 OXX=2412,9242,2445,4545",
			"GTAGCTTTCGACCGCTAGGCTAGCCCGAGATCGC"
		);

        REQUIRE(verify_id_wrapper(seq.name, 0).compare("2412") == 0);
		REQUIRE(verify_id_wrapper(seq.name, 0).compare("2412") == 0);
		REQUIRE(verify_id_wrapper(seq.name, 1).compare("9242") == 0);
		REQUIRE(verify_id_wrapper(seq.name, 2).compare("2445") == 0);
		REQUIRE(verify_id_wrapper(seq.name, 3).compare("4545") == 0);
	}
	SECTION("Testing verification of ID type with map of species IDs")
	{
		sequence seq(
			">ID=HNAC_EYY AC=KANMEA_R45 OXX=2412,9242,2445,4545",
			"GTAGCTTTCGACCGCTAGGCTAGCCCGAGATCGC"
		);

		std::unordered_map<std::string, std::string> id_map = {
			{">ID=HNAC_EYY AC=KANMEA_R45 OXX=2412,9242,2445,4545", "8892"},
            {">ID=YACD_EYY AC=NALMEA_R46 OXX=2563,9281,5823,4829", "7183"},
            {">ID=AHHC_IWO AC=YYYTWA_B25 OXX=0131,0992,8819,1121", "0914"}
		};

		REQUIRE(link_mod.verify_id_type(seq.name, &id_map).compare("8892") == 0);
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
    /* TODO: fix error in integration section
    SECTION( "Verifying metadata file usage feature successfully performs and integrates with link module." )
    {
        std::cout << "\n\n\nBeginning of metadata file usage test...\n";
        module_link mod;
        options_link opts;
        opts.metadata_fname =  "../test/input_data/full_design_clean_min30_taxtweak_100perc_jingmens_2019-09-12_segment.metadata,Name,Species";
        opts.prot_file_fname = "../test/input_data/full_design_clean_min30_taxtweak_100perc_jingmens_2019-09-12_segment.fasta";
        opts.peptide_file_fname = "../test/input_data/PV1_10K3000_53_encoded_segment.faa";
        opts.kmer_size = 7;
        opts.output_fname = "../test/test_metadata_output.txt";
        mod.run( &opts );
        REQUIRE( !opts.output_fname.empty() );
        std::cout << "End of metadata file usage test!\n\n\n";
    }
    */
}

TEST_CASE("Test zscore calculation", "[module_zscore]")
{
	SECTION("Run symmetrical trim option")
	{
		module_zscore z_mod;
		options_zscore z_opts;
		std::vector<nan_report> nan_values;
		peptide_score_data_sample_major input_data;

		z_opts.in_fname = "../test/input_data/test_zscore_score_matrix.tsv";
		z_opts.in_bins_fname = "../test/input_data/test_zscore_bins.tsv";
		z_opts.out_fname = "../test/test_zscore_trim_output.tsv";
		z_opts.hpd_percent = 0.00;
		z_opts.trim_percent = 0.90;
		
		std::ifstream bins_file(z_opts.in_bins_fname, std::ios_base::in);

		peptide_scoring::parse_peptide_scores(input_data, z_opts.in_fname);
		bin_collection peptide_bins = peptide_bin_io::parse_bins(bins_file);

		auto& zscore_matrix = input_data.scores;
		for (const auto sample_pair : zscore_matrix.get_row_labels())
		{
			std::string sample_name = sample_pair.first;
			z_mod.calculate_zscores(
				peptide_bins, &z_opts, zscore_matrix,
				sample_name, std::back_inserter(nan_values)
			);
		}
		peptide_scoring::write_peptide_scores(z_opts.out_fname, input_data);
	}
	SECTION("Verify trim option properly trims both head and tail")
	{
		std::ifstream ifexpected(
			"../test/expected/test_expected_zscore_trim_output.tsv",
			std::ios_base::in
		);
		std::ifstream ifactual(
			"../test/test_zscore_trim_output.tsv",
			std::ios_base::in
		);
		std::string expected_line;
		std::string actual_line;

		while (!ifexpected.eof() && !ifactual.eof())
		{
			std::getline(ifexpected, expected_line);
			std::getline(ifactual, actual_line);
			REQUIRE(expected_line.compare(actual_line) == 0);
		}
	}
	SECTION("Run highest density interval (hdi) option")
	{
		module_zscore z_mod;
		options_zscore z_opts;
		std::vector<nan_report> nan_values;
		peptide_score_data_sample_major input_data;

		z_opts.in_fname = "../test/input_data/test_zscore_score_matrix.tsv";
		z_opts.in_bins_fname = "../test/input_data/test_zscore_bins.tsv";
		z_opts.out_fname = "../test/test_zscore_hdi_output.tsv";
		z_opts.hpd_percent = 0.75;
		z_opts.trim_percent = 0.00;
		
		std::ifstream bins_file(z_opts.in_bins_fname, std::ios_base::in);

		peptide_scoring::parse_peptide_scores(input_data, z_opts.in_fname);
		bin_collection peptide_bins = peptide_bin_io::parse_bins(bins_file);

		auto& zscore_matrix = input_data.scores;
		for (const auto sample_pair : zscore_matrix.get_row_labels())
		{
			std::string sample_name = sample_pair.first;
			z_mod.calculate_zscores(
				peptide_bins, &z_opts, zscore_matrix,
				sample_name, std::back_inserter(nan_values)
			);
		}
		peptide_scoring::write_peptide_scores(z_opts.out_fname, input_data);
	}
	SECTION("Verify high density interval option properly trimmed distribution")
	{
		std::ifstream ifexpected(
			"../test/expected/test_expected_zscore_hdi_output.tsv",
			std::ios_base::in
		);
		std::ifstream ifactual(
			"../test/test_zscore_hdi_output.tsv",
			std::ios_base::in
		);
		std::string expected_line;
		std::string actual_line;

		while (!ifexpected.eof() && !ifactual.eof())
		{
			std::getline(ifexpected, expected_line);
			std::getline(ifactual, actual_line);
			REQUIRE(expected_line.compare(actual_line) == 0);
		}
	}
}

