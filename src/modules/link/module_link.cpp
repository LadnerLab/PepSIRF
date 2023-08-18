#include "logger.h"
#include "module_link.h"
#include "kmer_tools.h"
#include "fasta_parser.h"
#include "omp_opt.h"

#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/regex.hpp>
#include <boost/algorithm/string/join.hpp>
#include <boost/range/adaptor/transformed.hpp>


module_link::module_link()
{
    name = "Link";
}

void module_link::run( options *opts )
{
    options_link *l_opts = (options_link*) opts;

    fasta_parser fp;

    std::vector<sequence> peptides = fp.parse( l_opts->peptide_file_fname );
    std::vector<sequence> proteins = fp.parse( l_opts->prot_file_fname    );

    std::unordered_map<
        std::string,
        std::unordered_set<scored_entity<std::string,double>>
    > kmer_sp_map;

    std::vector<
        std::tuple<
            std::string,std::unordered_set<
                scored_entity<std::string,double>
            >
        >
    > peptide_sp_map;
    // When a metadata file is provided it is parsed for specific data and a
    // kmer map is created, otherwise process uses ID index to create a kmer
    // map
    if( l_opts->metadata_fname.empty() )
        {
            create_prot_map( kmer_sp_map, proteins, l_opts->k, l_opts->id_index );
        }
    else
        {
            Log::warn(
                "Metadata file has been provided and will be implemented with"
                " taxonomic ID index ignored\n"
            );
            metadata_map mp = metadata_map();
            std::unordered_map<std::string, std::string> meta_map;
            mp.build_map( &meta_map, l_opts->metadata_fname );
            create_prot_map( kmer_sp_map, proteins, l_opts->k, &meta_map );
        }

    if( l_opts->penalize_kmers )
        {
            create_pep_map_with_kmer_penalty( kmer_sp_map,
                                              peptide_sp_map,
                                              peptides,
                                              l_opts->k
                                              );
        }
    else
        {
            create_pep_map( kmer_sp_map,
                            peptide_sp_map,
                            peptides,
                            l_opts->k
                            );
        }

    write_outputs( l_opts->output_fname, peptide_sp_map );
}

void module_link::write_outputs(
    std::string fname,
    std::vector<
        std::tuple<
            std::string,
            std::unordered_set<scored_entity<std::string,double>>
        >
    >& peptide_sp_vec
) {
    std::ofstream out_file( fname );

    out_file << "Peptide Name\tLinked Species IDs with counts\n";

    for( auto it = peptide_sp_vec.begin(); it != peptide_sp_vec.end(); ++it ) // iterate over in vector
        {
            out_file << std::get<0>( *it ) << "\t"; // access tuple space 1

            auto in_vec = std::get<1>( *it ); // access tuple space 2

            std::size_t in_index = 0;
            for( auto it = in_vec.begin(); it != in_vec.end(); ++it )
                {
                        out_file << it->get_key()
                                 << ":" << it->get_score();

                    if( in_index != in_vec.size() - 1 )
                        {
                            out_file << ",";
                        }
                    ++in_index;
                }
            out_file << "\n";
        }
}

template<typename retrieve_id>
void module_link::create_prot_map(
    std::unordered_map<
        std::string,
        std::unordered_set<scored_entity<std::string,double>>
    >& scores_map,
    std::vector<sequence>& sequences,
    std::size_t k,
    retrieve_id retriever
) {
    std::size_t index   = 0;
    std::string spec_id = "";

    double t_start = omp_get_wtime();
    size_t num_prot = 0;

    std::size_t missing_spec_id_count = 0;
    std::ofstream ex_seqs_log(
        "excluded_protein_sequences.txt", std::ios::out
    );

    for( index = 0; index < sequences.size(); ++index )
        {
            ++num_prot;
            std::vector<std::string> kmers;
            kmer_tools::get_kmers( kmers, sequences[ index ].seq, k );
            std::unordered_set<scored_entity<std::string,double>>
                val_set;

            spec_id = verify_id_type( sequences[ index ].name, retriever );

            // check value for species ID column was missing in metadata file
            if (spec_id.empty())
            {
                ex_seqs_log << sequences[index].name << "\n";
                missing_spec_id_count += 1;

                kmers.clear();
                continue; // consider next sequence
            }

            for (auto it = kmers.begin(); it != kmers.end(); ++it)
            {
                // only inserts if key not already in map
                auto pair = std::get<0>(
                    scores_map.insert( std::make_pair( *it, val_set ) )
                );

                auto scored_ent = pair->second.emplace(
                    scored_entity<std::string,double>( spec_id, 0 )
                ).first;
                ++( scored_ent->get_score() );
            }

            // TODO: remove? since the vector is realloced every iteration
            // (line 130)
            kmers.clear();
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
    return metadata_map::get_id(sequence_data, map);
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

