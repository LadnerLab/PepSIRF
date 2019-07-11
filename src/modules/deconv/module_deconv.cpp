#include <iostream>
#include <fstream>
#include <sstream>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

#include "module_deconv.h"

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

    auto pep_species_vec = parse_linked_file( d_opts->linked_fname );


    sequential_map<std::size_t, std::vector<std::string>> id_pep_map;
    sequential_map<std::string, std::vector<std::size_t>> pep_id_map;
    sequential_map<std::size_t, std::size_t> id_count_map;

    id_to_pep( id_pep_map, pep_species_vec );
    pep_to_id( pep_id_map, pep_species_vec );
    count_species( id_count_map, pep_species_vec );

    std::cout << "Map size " << pep_species_vec.size() << "\n";

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

void module_deconv::count_species( sequential_map<std::size_t, std::size_t>&
                                   id_count_map,
                                   std::vector<std::pair<std::string, std::vector<std::size_t>>>&
                                   vector
                                 )
{
    for( auto it = vector.begin(); it != vector.end(); ++it )
        {
            std::vector<std::size_t>& id_ref = std::get<1>( *it );

            for( auto in = id_ref.begin(); in != id_ref.end(); ++in )
                {
                    auto count = id_count_map.find( *in );
                    if( count == id_count_map.end() )
                        {
                            id_count_map.emplace( *in, 0 );
                        }
                    id_count_map[ *in ]++;
                }
        }
}
