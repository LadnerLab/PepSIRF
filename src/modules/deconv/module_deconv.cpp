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

    // auto id_pep_map = 


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

                            // auto spec_vec
                            //     = std::make_shared<std::vector<std::shared_ptr<species_data>>>
                            //     ( std::vector<std::shared_ptr<species_data>>() );
                            // for( auto it = comma_delimited.begin(); it != comma_delimited.end(); ++it )
                            //     {
                            //         std::size_t id = boost::lexical_cast<std::size_t>( *it );

                            //         if( ret_map.find( id ) == ret_map.end() ) 
                            //             {
                            //                 ret_map.insert( std::make_pair( id,
                            //                                                 std::make_shared<species_data>
                            //                                                 ( species_data() )
                            //                                               )
                            //                               );
                            //                 spec_vec->push_back( ret_map.find( id )->second );
                            //             }
                            //         auto item = ret_map.find( id );
                            //         item->second->count += get_score( comma_delimited.size() );
                            //         item->second->common_peptides = spec_vec;
                            //     }
                        }
                }
            ++lineno;
        }

    return ret_vec;
}

int module_deconv::get_score( std::size_t size )
{
    return 1;
}
