#include <iostream>
#include <fstream>
#include <sstream>
#include <boost/algorithm/string.hpp>
#include <omp.h>
#include <boost/lexical_cast.hpp>

#include "module_deconv.h"
#include "time_keep.h"

module_deconv::module_deconv() = default;

species_data::species_data() 
{
    count = 0;
}

std::string module_deconv::get_name()
{
    return "Deconv";
}

void module_deconv::run( options *opts )
{
    options_deconv *d_opts = ( options_deconv * ) opts;
    time_keep::timer timer;

    timer.start = omp_get_wtime();

    auto pep_species_vec = parse_linked_file( d_opts->linked_fname );
    std::size_t thresh = d_opts->threshold;

    omp_set_num_threads( 2 );


    sequential_map<std::size_t, std::vector<std::string>> id_pep_map;
    sequential_map<std::string, std::vector<std::size_t>> pep_id_map;

    // vector holding the species ids with the highest count
    std::vector<std::pair<std::size_t, std::size_t>> id_counts;


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


    count_species( id_counts, id_pep_map );
    filter_counts( id_counts, thresh );
    std::cout << "Species id\tCount\n";

    while( id_counts.size() )
        {
            pep_species_vec.clear();

            auto max = id_counts.begin();
            auto max_id = std::get<0>( *max );

            std::cout << max_id << "\t" << std::get<1>( *max ) << "\n";

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


            count_species( id_counts, id_pep_map );
            // recreate id_counts
            filter_counts( id_counts, thresh );
        }

    timer.end = omp_get_wtime();

    std::cout << "Took " << time_keep::get_elapsed( timer ) << " seconds.\n";



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

int module_deconv::get_score( std::size_t size )
{
    return 1;
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

void module_deconv::count_species( std::vector<std::pair<std::size_t, std::size_t>>&
                                   id_counts,
                                   sequential_map<std::size_t,std::vector<std::string>>&
                                   id_count_map
                                 )
{
    id_counts.reserve( id_count_map.size() );

    for( auto it = id_count_map.begin(); it != id_count_map.end(); ++it )
        {
            id_counts.emplace_back( it->first, it->second.size() );
        }
    std::sort( id_counts.begin(), id_counts.end(),
               compare_pair<std::size_t>()
             );
}

void module_deconv::filter_counts( std::vector<std::pair<std::size_t, std::size_t>>& id_counts,
                                   std::size_t thresh
                                 )
{
    auto it = std::remove_if( id_counts.begin(), id_counts.end(),
                              [&]( std::pair<std::size_t, std::size_t> i )
                              { return std::get<1>( i ) < thresh; }
                            );
    id_counts.erase( it, id_counts.end() );
}
