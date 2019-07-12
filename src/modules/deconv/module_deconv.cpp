#include <iostream>
#include <fstream>
#include <sstream>
#include <boost/algorithm/string.hpp>
#include <omp.h>
#include <boost/lexical_cast.hpp>

#include "module_deconv.h"
#include "time_keep.h"

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

    auto enriched_species = parse_enriched_file( d_opts->enriched_fname );
    auto pep_species_vec  = parse_linked_file( d_opts->linked_fname );

    // filter out the peptides that are not enriched
    auto it = std::remove_if( pep_species_vec.begin(), pep_species_vec.end(),
                              [&]( std::pair<std::string,std::vector<std::size_t>>& i) -> bool 
                              { return enriched_species.find( std::get<0>( i ) )
                                      == enriched_species.end();
                              }
                            );
    pep_species_vec.erase( it, pep_species_vec.end() );

    double thresh = d_opts->threshold;

    omp_set_num_threads( d_opts->single_threaded ? 1 : 2 );


    sequential_map<std::size_t, std::vector<std::string>> id_pep_map;
    sequential_map<std::string, std::vector<std::size_t>> pep_id_map;

    // vector holding the species ids with the highest count
    std::vector<std::pair<std::size_t, double>> id_counts;

    std::vector<std::pair<std::size_t, double>> output_counts;

    #pragma omp parallel sections
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


    count_species( id_counts, id_pep_map, pep_id_map,
                   score_method::score_strategy::INTEGER_SCORING
                 );
    filter_counts( id_counts, thresh );

    while( id_counts.size() )
        {
            pep_species_vec.clear();

            auto max = id_counts.begin();
            double max_id = std::get<0>( *max );

            output_counts.emplace_back( max_id, std::get<1>( *max ) );

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
            id_counts.clear();


            #pragma omp parallel sections
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


            count_species( id_counts, id_pep_map, pep_id_map,
                           score_method::score_strategy::INTEGER_SCORING
                         );

            // recreate id_counts
            filter_counts( id_counts, thresh );
        }

    timer.end = omp_get_wtime();

    std::cout << "Took " << time_keep::get_elapsed( timer ) << " seconds.\n";

    write_outputs( d_opts->output_fname, output_counts );



}

std::vector<std::pair<std::string,std::vector<std::size_t>>>
module_deconv::parse_linked_file( std::string fname )
{
    std::ifstream in_stream( fname );

    
    std::vector<std::pair<std::string,std::vector<std::size_t>>>
        ret_vec;

    std::string line;

    int lineno = 0;

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

                            std::vector<std::size_t> id_ints;
                            std::for_each( comma_delimited.begin(), comma_delimited.end(),
                                           [&]( const std::string& item )
                                           { id_ints.push_back(
                                             boost::lexical_cast<std::size_t>( item )
                                                               );
                                           }
                                         );
                            ret_vec.emplace_back( split_line[ 0 ], id_ints );

                        }
                }
            ++lineno;
        }

    return ret_vec;
}

void module_deconv::id_to_pep( sequential_map<std::size_t, std::vector<std::string>>&
                               id_pep_map,
                               std::vector<std::pair<std::string, std::vector<std::size_t>>>&
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
                    auto find = id_pep_map.find( *inner );
                    if( find == id_pep_map.end() )
                        {
                            find = std::get<0>(
                                   id_pep_map.emplace( *inner, std::vector<std::string>() )
                                              );
                        }
                    find->second.push_back( pep );
                }
        }
}

double module_deconv::get_score( sequential_map<std::string,std::vector<std::size_t>>&
                                 spec_count_map,
                                 std::vector<std::string>& peptides,
                                 score_method::score_strategy strat
                               )
{
    if( strat != score_method::score_strategy::INTEGER_SCORING )
        {
            return (double) peptides.size();
        }

    double score = 0;

    std::size_t index = 0;

    #pragma omp parallel for private( index ) \
            shared( spec_count_map ) reduction( +:score )
    for( index = 0; index < peptides.size(); ++index )
        {
            score += 1.0 / (double) spec_count_map[ peptides[ index ] ].size();
        }
    return score;
}
void module_deconv::pep_to_id( sequential_map<std::string, std::vector<std::size_t>>&
                               pep_id_map,
                               std::vector<std::pair<std::string, std::vector<std::size_t>>>&
                               pep_species_vec
                             )
{
    for( auto it = pep_species_vec.begin(); it != pep_species_vec.end(); ++it )
        {
            pep_id_map.emplace( std::get<0>( *it ), std::get<1>( *it ) );
        }
}

void module_deconv::count_species( std::vector<std::pair<std::size_t, double>>&
                                   id_counts,
                                   sequential_map<std::size_t,std::vector<std::string>>&
                                   id_count_map,
                                   sequential_map<std::string,std::vector<std::size_t>>&
                                   spec_count_map,
                                   score_method::score_strategy strat
                                 )
{
    id_counts.reserve( id_count_map.size() );

    for( auto it = id_count_map.begin(); it != id_count_map.end(); ++it )
        {
            double score = get_score( spec_count_map,
                                      it->second,
                                      strat
                                    );
            id_counts.emplace_back( it->first, score );
        }
    std::sort( id_counts.begin(), id_counts.end(),
               compare_pair<std::size_t, double>()
             );

}

void module_deconv::filter_counts( std::vector<std::pair<std::size_t, double>>& id_counts,
                                   double thresh
                                 )
{
    auto it = std::remove_if( id_counts.begin(), id_counts.end(),
                              [&]( std::pair<std::size_t, double> i )
                              {
                                  return std::get<1>( i ) < thresh;
                              } 
                            );
    id_counts.erase( it, id_counts.end() );
}

void module_deconv::write_outputs( std::string out_name,
                                   std::vector<std::pair<std::size_t,double>>& out_counts
                                 )
{
    std::ofstream out_file( out_name );

    out_file << "Species ID\tCount\n";

    std::for_each( out_counts.begin(), out_counts.end(),
                   [&]( std::pair<std::size_t,double> elem )
                   {
                       out_file << std::get<0>( elem ) << "\t" <<
                           std::get<1>( elem ) << "\n";
                   }
                 );

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
