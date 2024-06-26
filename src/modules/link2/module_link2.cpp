#include "logger.h"
#include "module_link2.h"
#include "kmer_tools.h"
#include "fasta_parser.h"
#include "omp_opt.h"

#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/regex.hpp>
#include <boost/algorithm/string/join.hpp>
#include <boost/range/adaptor/transformed.hpp>

module_link2::module_link2()
{
    name = "Link2";
}

void module_link2::run( options *opts )
{
    struct time_keep::timer timer;
    timer.start();

    options_link2 *l_opts = (options_link2*) opts;

    fasta_parser fp;
    metadata_map2 mp;

    // read in metadata
    std::unordered_map<std::string, std::unordered_set<std::string>> meta_map;
    mp.build_map( &meta_map, l_opts->metadata_info );

    // read in peptide sequences
    std::vector<sequence> peptides = fp.parse( l_opts->peptide_file_fname );
    // read in target sequences
    std::vector<sequence> proteins = fp.parse( l_opts->prot_file_fname );

    std::vector<std::vector<std::size_t>> all_patts = generate_patterns(l_opts->kmer_size, l_opts->span);

    // TODO: implement parallelization
    omp_set_num_threads( l_opts->num_threads );

    std::unordered_map<std::string, std::unordered_set<std::string>>::iterator meta_idx;
    #pragma omp parallel for private( meta_idx ) shared( peptides, proteins, meta_map )

    #pragma omp parallel shared(peptides, proteins, meta_map)
    {
        std::unordered_map<std::string, std::unordered_map<std::string, int>> out_scores;

        #pragma omp for
        for( meta_idx = meta_map.begin(); meta_idx != meta_map.end(); meta_idx++ )
        {
            std::unordered_map<std::string, int> peptide_scores;

            // Log::info(meta_idx->first + "\n");

            for( const auto &patt : all_patts)
            {
                /*
                Log::info("Pattern: ");
                std::size_t temp_ct = 0;
                for( const auto& val : patt )
                {
                    Log::info(std::to_string(val));
                    if( temp_ct < patt.size() - 1 )
                    {
                        Log::info(", ");
                    }
                    temp_ct++;
                }
                Log::info("\n");
                */

                // generate all target kmers (for names in the metadata that are in the protein file) for this pattern
                std::unordered_set<std::string> target_kmers;
                std::unordered_set<std::string>::iterator meta_name;
                for( meta_name = (meta_idx->second).begin(); meta_name != (meta_idx->second).end(); meta_name++ )
                {
                    std::vector<sequence>::iterator protein_it;
                    protein_it = std::find_if( proteins.begin(), proteins.end(), 
                        [&](const sequence &s) {
                            return s.name == *meta_name;
                        } );
                    if( protein_it != proteins.end() )
                    {
                        get_patterned_kmers( target_kmers, protein_it->seq, patt, FILTER );
                    }
                }

                std::vector<sequence>::iterator peptide_it;
                for( peptide_it = peptides.begin(); peptide_it != peptides.end(); peptide_it++ )
                {
                    std::unordered_set<std::string> pep_kmers;
                    std::unordered_set<std::string> inter_kmers;
                    int score = 0;
                    get_patterned_kmers( pep_kmers, peptide_it->seq, patt, FILTER );

                    for( const auto& kmer : pep_kmers)
                    {
                        if( target_kmers.find(kmer) != target_kmers.end() )
                        {
                            score++;
                        }
                    }

                    // update score if found higher than previous pattern
                    // Note: accessing an undefined key in map will set it to 0
                    if( score > peptide_scores[ peptide_it->name ] )
                    {
                        peptide_scores[ peptide_it->name ] = score;
                    }
                }
            }

            // add the peptide scores to output map
            for( const auto& pep_score : peptide_scores )
            {
                if( pep_score.second > 0 )
                {
                    #pragma omp critical
                    {
                        out_scores[ pep_score.first ][ meta_idx->first ] = pep_score.second;
                    }
                }
            }
        }

        write_outputs( l_opts->output_fname, out_scores );
    }

    timer.stop();
    Log::info(
        "Took " + std::to_string(time_keep::get_elapsed(timer))
        + " second(s).\n"
    );
}

void module_link2::write_outputs(
    std::string fname,
    const std::unordered_map<std::string, std::unordered_map<std::string, int>> out_scores
) {
    std::ofstream out_file( fname );

    out_file << "Peptide Name\tLinked Species IDs with counts\n";

    for( auto it = out_scores.begin(); it != out_scores.end(); ++it ) // iterate over in vector
        {
            out_file << it->first << "\t"; // access peptide name

            auto in_vec = it->second; // access score vector

            std::size_t in_index = 0;
            for( auto it = in_vec.begin(); it != in_vec.end(); ++it )
                {
                        out_file << it->first
                                 << ":" << it->second;

                    if( in_index != in_vec.size() - 1 )
                        {
                            out_file << ",";
                        }
                    ++in_index;
                }
            out_file << "\n";
        }
}

void module_link2::get_patterned_kmers(std::unordered_set<std::string> &target_kmers, std::string seq, 
                                        std::vector<std::size_t> patt, std::vector<char> filter)
{
    std::size_t span = patt.size();
    std::size_t seq_idx;
    std::size_t patt_idx;
    bool filter_out = false;

    for( seq_idx = 0; seq_idx < (seq.size() - span + 1); seq_idx++ )
    {
        std::string curr_kmer = seq.substr(seq_idx, span);
        std::string patt_kmer;
        for( const auto& filt_let : filter )
        {   
            if( curr_kmer.find(filt_let) != std::string::npos )
            {
                filter_out = true;
            }
        }
        if( !filter_out )
        {
            for( patt_idx = 0; patt_idx < patt.size(); patt_idx++ )
            {
                if( patt[patt_idx] )
                {
                    patt_kmer.push_back(curr_kmer[patt_idx]);
                }
            }
        }
        filter_out = false;

        target_kmers.insert(patt_kmer);
    }

}

std::vector<std::vector<std::size_t>> module_link2::generate_patterns( std::size_t kmer_size, std::size_t span )
{
    std::vector<std::vector<std::size_t>> all_patts;
    std::size_t iter;

    if( kmer_size > span )
    {
        Log::error("Kmer size is greater than span!");
    }

    std::vector<std::size_t> patt( span - kmer_size, 0 );
    std::vector<std::size_t> kmer_patt( kmer_size, 1 );
    patt.insert( patt.end(), kmer_patt.begin(), kmer_patt.end() );

    // add all permutations
    do {
        all_patts.emplace_back( patt );
    } while( std::next_permutation(patt.begin(), patt.end()) );

    return all_patts;
}