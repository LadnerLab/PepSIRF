#include <iostream>
#include <fstream>
#include <sstream>
#include <boost/algorithm/string.hpp>
#include <omp.h>
#include <boost/lexical_cast.hpp>
#include <boost/regex.hpp>

#include "module_deconv.h"
#include "kmer_tools.h"
#include "time_keep.h"
#include "fasta_parser.h"

module_deconv::module_deconv() = default;

std::string module_deconv::get_name()
{
    return "Deconv";
}

void module_deconv::run( options *opts )
{
    options_deconv *d_opts = ( options_deconv * ) opts;
    time_keep::timer timer;

    timer.start = omp_get_wtime();

    if( d_opts->create_linkage )
        {
            create_linkage( d_opts );
        }
    else 
        {
            choose_kmers( d_opts );
        }

    timer.end = omp_get_wtime();

    std::cout << "Took " << time_keep::get_elapsed( timer ) << " second(s).\n";
}

void module_deconv::choose_kmers( options_deconv *opts )
{
    options_deconv *d_opts = opts;

    auto enriched_species = parse_enriched_file( d_opts->enriched_fname );
    auto pep_species_vec  = parse_linked_file( d_opts->linked_fname );

    score_method::score_strategy score_strat   = get_score_method( d_opts );
    score_method::filter_strategy filter_strat = get_filter_method( d_opts );

    // filter out the peptides that are not enriched
    auto it = std::remove_if( pep_species_vec.begin(), pep_species_vec.end(),
                              [&]( std::pair<std::string,std::vector<std::pair<std::size_t,std::size_t>>>& i) -> bool 
                              { return enriched_species.find( std::get<0>( i ) )
                                      == enriched_species.end();
                              }
                            );
    pep_species_vec.erase( it, pep_species_vec.end() );

    std::size_t thresh = d_opts->threshold;

    omp_set_num_threads( d_opts->single_threaded ? 1 : 2 );

    sequential_map<std::size_t, std::vector<std::string>> id_pep_map;
    sequential_map<std::string, std::vector<std::pair<std::size_t,std::size_t>>>
        pep_id_map;

    // vector holding the species ids with the highest count
    std::vector<std::pair<std::size_t, double>> species_scores;
    sequential_map<std::size_t, std::size_t> species_peptide_counts;

    std::vector<std::tuple<std::size_t, std::size_t, double>> output_counts;

    #pragma omp parallel
    {
        #pragma omp sections
        {
            #pragma omp section
            {
                id_to_pep( id_pep_map, pep_species_vec );
            }

            #pragma omp section
            {
                pep_to_id( pep_id_map, pep_species_vec );
            }
        }

        // implicit barrier
        #pragma omp sections
        {
            #pragma omp section
            {
                score_species( species_scores, id_pep_map, pep_id_map,
                               score_strat
                             );
            }

            #pragma omp section
            {
                get_species_counts_per_peptide( id_pep_map,
                                                species_peptide_counts
                                                );

                // recreate species_scores
                filter_counts<std::size_t,std::size_t>
                    ( species_peptide_counts, thresh );

            }
        }

    }

    while( species_peptide_counts.size() )
        {
            pep_species_vec.clear();

            auto max = species_scores.begin();
            std::size_t max_id = std::get<0>( *max );

            output_counts.emplace_back( std::make_tuple(
                                                        max_id,
                                                        species_peptide_counts
                                                        .find( max_id )->second, // count for speices
                                                        std::get<1>( *max )
                                                       )
                                      );

            auto max_peptides = id_pep_map.find( max_id )->second;

            for( auto pep_id = max_peptides.begin();
                 pep_id != max_peptides.end();
                 ++pep_id )
                {
                    pep_id_map.erase( *pep_id );
                }

            for( auto it = pep_id_map.begin(); it != pep_id_map.end(); ++it )
                {
                    pep_species_vec.emplace_back( it->first, it->second );
                }

            id_pep_map.clear();
            pep_id_map.clear();
            species_scores.clear();
            species_peptide_counts.clear();

            #pragma omp parallel
            {
                #pragma omp sections
                {
                    #pragma omp section
                    {
                        id_to_pep( id_pep_map, pep_species_vec );
                    }

                    #pragma omp section
                    {
                        pep_to_id( pep_id_map, pep_species_vec );
                    }
                }

                // implicit barrier
                #pragma omp sections
                {
                    #pragma omp section
                    {
                        score_species( species_scores, id_pep_map, pep_id_map,
                                       score_strat
                                     );
                    }

                    #pragma omp section
                    {
                        get_species_counts_per_peptide( id_pep_map,
                                                        species_peptide_counts
                                                        );

                        // recreate species_scores
                        filter_counts<std::size_t,std::size_t>
                            ( species_peptide_counts, thresh );

                    }
                }

            }

        }

    std::string id_name_map_fname;
    std::map<std::size_t,std::string> name_id_map;
    std::map<std::size_t,std::string>* name_id_map_ptr = nullptr;
    if( d_opts->id_name_map_fname.compare( "" ) )
        {
            parse_name_map( d_opts->id_name_map_fname, name_id_map );
            name_id_map_ptr = &name_id_map;
        }
    write_outputs( d_opts->output_fname, name_id_map_ptr, output_counts );

}

void module_deconv::create_linkage( options_deconv *opts )
{
    options_deconv *d_opts = opts;

    fasta_parser fp;
    std::vector<sequence> peptides = fp.parse( d_opts->peptide_file_fname );
    std::vector<sequence> proteins = fp.parse( d_opts->prot_file_fname    );

    sequential_map<std::string,
                   sequential_map<std::size_t,std::size_t>>
    kmer_sp_map;

    std::vector<std::tuple<std::string,sequential_map<std::size_t,std::size_t>>>
        peptide_sp_map;

    create_prot_map( kmer_sp_map, proteins, d_opts->k );
    create_pep_map( kmer_sp_map, peptide_sp_map, peptides, d_opts->k );
    write_outputs( d_opts->output_fname, peptide_sp_map );

}

std::vector<std::pair<std::string,std::vector<std::pair<std::size_t,std::size_t>>>>
module_deconv::parse_linked_file( std::string fname )
{
    std::ifstream in_stream( fname );

    std::vector<std::pair<std::string,std::vector<std::pair<std::size_t,std::size_t>>>>
        ret_vec;

    std::string line;

    int lineno = 0;
    const boost::regex id_count_re{ "([0-9]+):{0,1}([0-9]*)" };


    while( std::getline( in_stream, line ) )
        {
            if( lineno )
                {
                    std::vector<std::string> split_line;
                    boost::trim_right( line );
                    boost::split( split_line, line,
                                  boost::is_any_of( "\t" )
                                );


                    std::vector<std::string> comma_delimited;
                    if( split_line.size() > 1 )
                        {
                            boost::split( comma_delimited, split_line[ 1 ],
                                          boost::is_any_of( "," )
                                        );

                            std::vector<std::pair<std::size_t,std::size_t>> id_ints;

                            for( auto& item : comma_delimited )
                                {
                                    boost::smatch match;
                                    boost::regex_search( item, match, id_count_re );

                                    // matched 'id', no count found
                                    if( match[ 1 ] != ""
                                        && match[ 2 ] == "" )
                                        {
                                            id_ints.emplace_back(
                                                                 std::make_pair(
                                                                 boost::lexical_cast<std::size_t>
                                                                 ( match[ 1 ] ),
                                                                 boost::lexical_cast<std::size_t>
                                                                 ( 1 )
                                                                                )
                                                                );
                                        }
                                    // matched 'id:count'
                                    else if( match[ 1 ] != ""
                                             && match[ 2 ] != "" ) 
                                        {
                                            id_ints.emplace_back(
                                                                 std::make_pair(
                                                                 boost::lexical_cast<std::size_t>
                                                                 ( match[ 1 ] ),
                                                                 boost::lexical_cast<std::size_t>
                                                                 ( match[ 2 ] )
                                                                                )
                                                                );
                                        }
                                    else
                                        {
                                            std::cout << match.length() << "\n";
                                            std:: cout << item << "\n";
                                        }

                                }
                            ret_vec.emplace_back( split_line[ 0 ], id_ints );

                        }
                }
            ++lineno;
        }

    return ret_vec;
}

void module_deconv::id_to_pep( sequential_map<std::size_t, std::vector<std::string>>&
                               id_pep_map,
                               std::vector<std::pair<std::string, std::vector<std::pair<std::size_t,std::size_t>>>>&
                               pep_species_vec )
{
    for( auto it = pep_species_vec.begin(); it != pep_species_vec.end(); ++it )
        {
            std::string& pep = std::get<0>( *it );

            for( auto inner = std::get<1>( *it ).begin();
                     inner != std::get<1>( *it ).end();
                 ++inner
               )
                {
                    auto find = id_pep_map.find( std::get<0>( *inner ) );
                    if( find == id_pep_map.end() )
                        {
                            find = std::get<0>(
                                               id_pep_map.emplace( std::get<0>( *inner ), std::vector<std::string>() )
                                              );
                        }
                    find->second.push_back( pep );
                }
        }
}

double module_deconv::get_score( sequential_map<std::string,std::vector<std::pair<std::size_t,std::size_t>>>&
                                 spec_count_map,
                                 std::size_t id,
                                 std::vector<std::string>& peptides,
                                 score_method::score_strategy strat
                               )
{
    double score = 0;

    if( strat == score_method::score_strategy::INTEGER_SCORING )
        {
            return (double) peptides.size();
        }

    else if( strat == score_method::score_strategy::FRACTIONAL_SCORING )
        {
            std::size_t index = 0;

            #pragma omp parallel for private( index )           \
             shared( spec_count_map ) reduction( +:score )
            for( index = 0; index < peptides.size(); ++index )
                {
                    score += 1.0 / (double) spec_count_map[ peptides[ index ] ].size();
                }
        }

    else if( strat == score_method::score_strategy::SUMMATION_SCORING )
        {
            std::size_t index = 0;

            #pragma omp parallel for private( index )           \
             shared( spec_count_map ) reduction( +:score )
            for( index = 0; index < peptides.size(); ++index )
                {
                    for( auto& iter : spec_count_map[ peptides[ index ] ] )
                        {
                            if( std::get<0>( iter ) == id )
                                {
                                    score += std::get<1>( iter );
                                }
                        }
                }
        }
    
    return score;
}

void module_deconv::pep_to_id( sequential_map<std::string, std::vector<std::pair<std::size_t,std::size_t>>>&
                               pep_id_map,
                               std::vector<std::pair<std::string, std::vector<std::pair<std::size_t,std::size_t>>>>&
                               pep_species_vec
                             )
{
    for( auto it = pep_species_vec.begin(); it != pep_species_vec.end(); ++it )
        {
            pep_id_map.emplace( std::get<0>( *it ),
                                std::get<1>( *it )
                              );
        }
}

void module_deconv::score_species( std::vector<std::pair<std::size_t, double>>&
                                   id_counts,
                                   sequential_map<std::size_t,std::vector<std::string>>&
                                   id_count_map,
                                   sequential_map<std::string,std::vector<std::pair<std::size_t,std::size_t>>>&
                                   spec_count_map,
                                   score_method::score_strategy strat
                                )

{
    id_counts.reserve( id_count_map.size() );

    for( auto it = id_count_map.begin(); it != id_count_map.end(); ++it )
        {
            double score = get_score( spec_count_map,
                                      it->first,
                                      it->second,
                                      strat
                                    );
            id_counts.emplace_back( it->first, score );
        }
    std::sort( id_counts.begin(), id_counts.end(),
               compare_pair<std::size_t, double>()
             );

}

void module_deconv::write_outputs( std::string out_name,
                                   std::map<std::size_t,std::string>*
                                   id_name_map,
                                   std::vector<
                                   std::tuple<std::size_t,std::size_t,double>
                                   >&
                                   out_counts
                                 )

{
    std::ofstream out_file( out_name );

    if( id_name_map != nullptr )
        {
            out_file << "Species Name\t";
        }


    out_file << "Species ID\tCount\tScore\n";


    for( auto& elem : out_counts )
        {
            if( id_name_map != nullptr
                && id_name_map->find( std::get<0>( elem ) ) != id_name_map->end() )
                {
                    out_file << id_name_map->find( std::get<0>( elem ) )->second << "\t";
                }
            else if( id_name_map != nullptr
                     && id_name_map->find( std::get<0>( elem ) ) == id_name_map->end() )
                {
                    out_file << "\t";
                }

            out_file << std::get<0>( elem ) << "\t" <<
                std::get<1>( elem ) << "\t" << 
                std::get<2>( elem ) << "\n";
        }

}

sequential_set<std::string>
module_deconv::parse_enriched_file( std::string f_name )
{
    sequential_set<std::string>
        return_val;

    std::ifstream in_file( f_name );
    std::string line;

    while( std::getline( in_file, line ) )
        {
            boost::trim_right( line );
            return_val.insert( line );
        }
    return return_val;
}

std::size_t module_deconv::get_id( std::string name )
{
    static const boost::regex id_re{ "OXX=([0-9]+),([0-9]*),([0-9]*),([0-9])" };
    boost::smatch match;

    if( boost::regex_search( name, match, id_re ) )
        {
            return boost::lexical_cast<std::size_t>( match[ 2 ] );
        }
    return 0;
}

void module_deconv::create_prot_map( sequential_map<std::string,
                                     sequential_map<std::size_t,std::size_t>>&
                                    map,
                                     std::vector<sequence>& sequences,
                                     std::size_t k
                                   )
{
    std::size_t index   = 0;
    std::size_t spec_id = 0;

    double t_start = omp_get_wtime();
    double num_prot = 0;

    for( index = 0; index < sequences.size(); ++index )
        {
            ++num_prot;
            std::vector<std::string> kmers;

            spec_id = get_id( sequences[ index ].name );
            kmer_tools::get_kmers( kmers, sequences[ index ].seq, k );
            sequential_map<std::size_t,std::size_t> val_map;

            for( auto it = kmers.begin(); it != kmers.end(); ++it )
                {
                    // only inserts if key not already in map
                    auto pair = std::get<0>(
                                            map.insert( std::make_pair( *it,
                                                                        val_map
                                                                      )
                                                       )
                                            );
                    auto pair_s = std::get<0>(
                                              pair->second.insert( std::make_pair( spec_id, 0 ) )
                                             );
                    ++(pair_s->second);
                }

            kmers.clear();
        }

    double t_end = omp_get_wtime();

    std::cout << num_prot << " proteins done in " << t_end - t_start << " seconds. (" << ( t_end - t_start ) / num_prot << " seconds per peptide\n";
}

void module_deconv::create_pep_map( sequential_map<std::string,
                                    sequential_map<std::size_t,std::size_t>>&
                                    kmer_sp_map,
                                    std::vector<std::tuple<std::string,sequential_map<std::size_t,std::size_t>>>&
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
            sequential_map<std::size_t,std::size_t> ids;

            kmer_tools::get_kmers( kmers, peptides[ index ].seq, k );

            peptide_sp_vec.insert( peptide_sp_vec.begin() + index, std::make_tuple( peptides[ index ].name, ids ) );

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
                                    if( id_ref.find( it->first ) == id_ref.end() )
                                        {
                                            id_ref.insert( std::make_pair( it->first, 1 ) );
                                        }
                                    else
                                        {
                                            id_ref[ it->first ] += 1;
                                        }
                                }
                        }
                }
            kmers.clear();
        }
}

void module_deconv::write_outputs( std::string fname,
                                   std::vector<std::tuple<std::string,sequential_map<std::size_t,std::size_t>>>&
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
                    out_file << std::get<0>( *it )
                             << ":" << std::get<1>( *it );

                    if( in_index != in_vec.size() - 1 )
                        {
                            out_file << ",";
                        }
                    ++in_index;
                }
            out_file << "\n";
        }
}


score_method::filter_strategy
module_deconv::get_filter_method( options_deconv *opts )
{
    if( opts->count_filtering )
        {
            return score_method::filter_strategy::COUNT_FILTER;
        }
    return score_method::filter_strategy::SCORE_FILTER;
}

score_method::score_strategy
module_deconv::get_score_method( options_deconv *opts )
{
    if( !( opts->fractional_scoring
           || opts->summation_scoring
          )
      )
        {
            return score_method::score_strategy::INTEGER_SCORING;
        }

    return  opts->fractional_scoring ?
             score_method::score_strategy::FRACTIONAL_SCORING :
             score_method::score_strategy::SUMMATION_SCORING;
}


void
module_deconv::parse_name_map( std::string fname,
                               std::map<std::size_t,std::string>& name_map
                             )
{
    std::ifstream in_stream( fname );
    std::size_t lineno = 0;

    std::string line;

    while( std::getline( in_stream, line ) )
        {
            std::vector<std::string> split_line;
            if( lineno )
                {

                    boost::trim_right( line );
                    boost::split( split_line, line,
                                  boost::is_any_of( "|" )
                                );

                    // if this is a species-level ID,
                    // then this will provide the mapping.
                    // Otherwise some space is wasted, but that's not
                    // a major concern.
                    boost::trim( split_line[ 0 ] );
                    boost::trim( split_line[ 1 ] );

                    std::size_t id = boost::lexical_cast<std::size_t>
                        ( split_line[ 0 ] );
                    std::string second = split_line[ 1 ];
                    name_map.insert( std::make_pair( id, second ) );
                }

            ++lineno;
        }
}


void module_deconv::
get_species_counts_per_peptide( sequential_map<std::size_t, std::vector<std::string>>&
                                id_pep_map,
                                sequential_map<std::size_t,std::size_t>& pep_counts
                              )
{
    for( auto count = id_pep_map.begin();
         count != id_pep_map.end();
         ++count
       )
        {
            pep_counts.emplace( count->first,
                                count->second.size()
                              );
        }
}
