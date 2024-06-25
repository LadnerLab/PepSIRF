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

    std::unordered_map<std::string, std::unordered_map<std::string, int>> out_scores;

    // TODO: implement parallelization
    std::unordered_map<std::string, std::unordered_set<std::string>>::iterator meta_idx;
    for( meta_idx = meta_map.begin(); meta_idx != meta_map.end(); meta_idx++ )
    {
        std::unordered_map<std::string, int> peptide_scores;

        // Log::info(meta_idx->first + "\n");

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
                out_scores[ pep_score.first ][ meta_idx->first ] = pep_score.second;
            }
        }
    }

    write_outputs( l_opts->output_fname, out_scores );

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

/*
template<typename retrieve_id>
void module_link::filter_proteins(
    std::vector<sequence> *proteins, 
    const std::unordered_map<std::string, std::string> meta_map
    ) {
}

    // check for at least one species ID missing for proteins
    if (missing_spec_id_count > 0)
    {
        std::cout << "WARNING: " << missing_spec_id_count
            << " sequence(s) in metadata file did not have a value for their"
            << " \"SpeciesID\" column - they have been excluded from this run.\n"
            << "Please review the aforementioned sequences in"
            << " \"excluded_protein_sequences.txt\"\n";
    }
    // otherwise, assume all specified sequences had a species ID
    else
    {
        // temp namespace
        namespace fs = boost::filesystem;

        // remove excluded protein seqs list file - not needed
        fs::path ex_file_path = "./excluded_protein_sequences.txt";

        // result not captured because file is guaranteed to exist; maybe this
        // can be refactored in such a way the file is not created unless
        // written to
        fs::remove(ex_file_path);
    }

    double t_end = omp_get_wtime();
    std::string str_interval = std::to_string(t_end - t_start);
    std::string str_rate = std::to_string((t_end - t_start) / (double)num_prot);

    Log::info(
        std::to_string(num_prot) + " proteins done in " + str_interval
        + " seconds." + " (" + str_rate + " seconds per peptide)\n"
    );
}

std::string module_link::get_id( std::string name, std::size_t id_index )
{
    static const boost::regex id_re{ "OXX=(\\S*),(\\S*),(\\S*),(\\S*)" };

    boost::smatch match;

    if( boost::regex_search( name, match, id_re ) )
        {
            return match[ id_index + 1 ];
        }
    return 0;
}

std::string module_link::verify_id_type(std::string sequence_data, std::size_t id_index)
{
    return get_id(sequence_data, id_index);
}

std::string module_link::verify_id_type(std::string sequence_data, std::unordered_map<std::string, std::string> *map)
{
    std::string found_id = metadata_map::get_id(sequence_data, map);

    if (found_id.empty())
    {
        Log::error(
            "Protein file contains sequences not represented in metadata file!"
        );
    }

    return found_id;
}

void module_link::create_pep_map( std::unordered_map<std::string,
                                  std::unordered_set<scored_entity<std::string,double>>>&
                                  kmer_sp_map,
                                  std::vector<std::tuple<std::string,
                                  std::unordered_set<scored_entity<std::string,double>>>>&
                                  peptide_sp_vec,
                                  std::vector<sequence>&
                                  peptides,
                                  std::size_t k
                                )

{
    peptide_sp_vec.reserve( peptides.size() );

    std::size_t index = 0;
    for( index = 0; index < peptides.size(); ++index )
        {
            // get the kmers from this peptide
            std::vector<std::string> kmers;
            std::unordered_set<scored_entity<std::string,double>> ids;

            kmer_tools::get_kmers( kmers, peptides[ index ].seq, k );

            peptide_sp_vec.insert( peptide_sp_vec.begin() + index,
                                   std::make_tuple( peptides[ index ].name, ids )
                                 );

            std::size_t kmer_index = 0;

            auto& id_ref =
                std::get<1>( peptide_sp_vec[ index ] );

            // for each of this peptide's kmers grab the counts from kmer_sp_map
            for( kmer_index = 0; kmer_index < kmers.size(); ++kmer_index )
                {
                    auto id_count = kmer_sp_map.find( kmers[ kmer_index ] );
                    if( id_count != kmer_sp_map.end() )
                        {
                            for( auto it = id_count->second.begin(); it != id_count->second.end(); ++it )
                                {
                                    if( id_ref.find( *it ) == id_ref.end() )
                                        {
                                            id_ref.insert( scored_entity<std::string,double>
                                                           ( it->get_key(), 1 )
                                                         );
                                        }
                                    else
                                        {
                                            id_ref.find( *it )->get_score()++;
                                        }
                                }
                        }
                }
            kmers.clear();
        }
}

void module_link
::create_pep_map_with_kmer_penalty( std::unordered_map<std::string,
                                    std::unordered_set<scored_entity<std::string,double>>>&
                                    kmer_sp_map,
                                    std::vector<std::tuple<std::string,
                                    std::unordered_set<scored_entity<std::string,double>>>>&
                                    peptide_sp_vec,
                                    std::vector<sequence>&
                                    peptides,
                                    std::size_t k
                                  )

{
    peptide_sp_vec.reserve( peptides.size() );

    std::unordered_set<scored_entity<std::string,std::size_t>> kmer_pep_frequency;

    kmer_tools::get_kmer_frequencies( kmer_pep_frequency, peptides, k );

    auto normalize_score = []( const std::size_t score,
                               const std::size_t penalty
                             )
        -> double
        {
            return ( (double) score ) / ( (double) (penalty ) );
        };

    std::size_t index = 0;
    for( index = 0; index < peptides.size(); ++index )
        {
            // get the kmers from this peptide
            std::vector<std::string> kmers;
            std::unordered_set<scored_entity<std::string,double>> ids;

            kmer_tools::get_kmers( kmers, peptides[ index ].seq, k );

            peptide_sp_vec.insert( peptide_sp_vec.begin() + index,
                                   std::make_tuple( peptides[ index ].name, ids )
                                 );

            std::size_t kmer_index = 0;

            auto& id_ref =
                std::get<1>( peptide_sp_vec[ index ] );

            // for each of this peptide's kmers grab the counts from kmer_sp_map
            for( kmer_index = 0; kmer_index < kmers.size(); ++kmer_index )
                {
                    auto id_count = kmer_sp_map.find( kmers[ kmer_index ] );
                    if( id_count != kmer_sp_map.end() )
                        {
                            scored_entity<std::string,std::size_t>
                                temp_score( kmers[ kmer_index ], 0 );

                            std::size_t score_penalty = kmer_pep_frequency
                                                          .find( temp_score )->get_score();

                            for( auto it = id_count->second.begin(); it != id_count->second.end(); ++it )
                                {
                                    if( id_ref.find( *it ) == id_ref.end() )
                                        {
                                            id_ref.insert( scored_entity<std::string,double>
                                                           ( it->get_key(),
                                                             normalize_score( 1, score_penalty )
                                                           )
                                                         );
                                        }
                                    else
                                        {
                                            const auto& score = *id_ref.find( *it );
                                            score.get_score() +=
                                                normalize_score( 1, score_penalty );
                                        }
                                }
                        }
                }
            kmers.clear();
        }
}
*/

