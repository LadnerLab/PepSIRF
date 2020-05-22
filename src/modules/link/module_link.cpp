#include "module_link.h"
#include "kmer_tools.h"
#include "fasta_parser.h"
#include "omp_opt.h"

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/regex.hpp>
#include <boost/algorithm/string/join.hpp>
#include <boost/range/adaptor/transformed.hpp>



module_link::module_link() = default;

void module_link::run( options *opts )
{
    options_link *l_opts = (options_link*) opts;

    fasta_parser fp;

    std::vector<sequence> peptides = fp.parse( l_opts->peptide_file_fname );
    std::vector<sequence> proteins = fp.parse( l_opts->prot_file_fname    );

    std::unordered_map<std::string,
                       std::unordered_set<scored_entity<std::string,double>>
                       >
        kmer_sp_map;

    std::vector<
        std::tuple<
            std::string,std::unordered_set<
                scored_entity<
                    std::string,double
                    >
                            >
            >
        >
        peptide_sp_map;
    // When a metadata file is provided it is parsed for specific data and a kmer map is created, otherwise process uses ID index to create a kmer map
    if( l_opts->metadata_fname.length() != 0 )
        {
        //is there a way to know if the index was provided or if it is the default? For now, it will be a general warning when metadata is provided.
            std::cout << "WARNING: Metadata file has been provided and will be implemented with taxonomic ID index ignored." << std::endl;
        /*
        parse metadata file to obtain column data, create protein linkage map
        metadata vector will contain three values, file name, one column name, a second column name

            STEP1: open metadata file, verify given column names are found in metadata file
            STEP2: loop through protein sequence vector
                for each protein sequence: -increment protein count
                                           -get list of kmers - k is length of kmer
                                           STEP3:
                                                for each kmer in list(kmer==substring):
                                                    -emplace a pair (kmer, value)
                                                    -set second, value, in pair to a scored entity with ?string value==second column name (eg species)? and double value 0
                                                    -increment double value once
                -clear kmer list
            Is this right? Right now this really just leaves out get_id - since we already have all the ids for the specified column
        */
            //find and verify existence and order for file name/path, pep seq name, spec name : eg. taxtweak_2019-09-12.metadata,Name,Species
            std::vector<std::string> metadata_options;
            metadata_options.reserve( 3 );
            boost::split( metadata_options, l_opts->metadata_fname, boost::is_any_of( "," ) );
            std::ifstream metadata_file( metadata_options[0], std::ios_base::in );
            if( !metadata_file.is_open() )
                {
                    throw std::runtime_error( "File could not be opened. Verify metadata file exists.\n" );
                }
            if( metadata_options.size() < 3 )
                {
                    throw std::runtime_error( "Missing required specifications for meta flag. Use \"--meta [file name],[sequence name],[taxonomic identity]\" format.\n");
                }
            std::string line;
            std::getline( metadata_file, line );
            std::vector<std::string>::iterator opt_iter;
            opt_iter = ++metadata_options.begin(); // step past file name
            std::for_each( opt_iter, metadata_options.end(), [ line ]( std::string header, std::size_t next_pos, std::size_t last_pos )
                    {
                        if( next_pos = line.find( header ) == std::string::npos )
                            {
                                throw std::runtime_error( "Header '" + header + "' could not be found in metadata file. Verify names being entered to names in metadata file.\n" );
                            }
                        if( last_pos > next_pos )
                            {
                                throw std::runtime_error( "Order of provided names do not match meta flag format. Use \"--meta [file name],[sequence name],[taxonomic identity]\" format.\n" );
                            }
                    }
                );
        }
    else
        {
            create_prot_map( kmer_sp_map, proteins, l_opts->k, l_opts->id_index );
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

void module_link::write_outputs( std::string fname,
                                 std::vector<std::tuple<std::string,
                                 std::unordered_set<scored_entity<std::string,double>>>>&
                                 peptide_sp_vec
                               )
{
    std::ofstream out_file( fname );

    out_file << "Peptide Name\tLinked Species IDs with counts\n";

    for( auto it = peptide_sp_vec.begin(); it != peptide_sp_vec.end(); ++it )
        {
            out_file << std::get<0>( *it ) << "\t";

            auto in_vec = std::get<1>( *it );

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

void module_link::create_prot_map( std::unordered_map<std::string,
                                   std::unordered_set<scored_entity<std::string,double>>>&
                                   scores_map,
                                   std::vector<sequence>& sequences,
                                   std::size_t k,
                                   std::size_t id_index
                                 )
{
    std::size_t index   = 0;
    std::string spec_id = "";

    double t_start = omp_get_wtime();
    double num_prot = 0;

    for( index = 0; index < sequences.size(); ++index )
        {
            ++num_prot;
            std::vector<std::string> kmers;

            spec_id = get_id( sequences[ index ].name, id_index );
            kmer_tools::get_kmers( kmers, sequences[ index ].seq, k );
            std::unordered_map<std::string,std::size_t> val_map;

            std::unordered_set<scored_entity<std::string,double>>
                val_set;

            for( auto it = kmers.begin(); it != kmers.end(); ++it )
                {
                    // only inserts if key not already in map
                    auto pair = std::get<0>(
                                            scores_map.insert( std::make_pair( *it,
                                                                               val_set
                                                                               )
                                                             )
                                            );

                    auto scored_ent = pair->second.emplace( scored_entity<std::string,double>
                                                            ( spec_id, 0 )
                                                          )
                                      .first;
                    ++( scored_ent->get_score() );
                }

            kmers.clear();
        }

    double t_end = omp_get_wtime();

    std::cout << num_prot << " proteins done in "
              << t_end - t_start
              << " seconds. (" << ( t_end - t_start ) / num_prot
              << " seconds per peptide\n";
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
