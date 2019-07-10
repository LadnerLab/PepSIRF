#include <iostream>
#include <fstream>
#include <sstream>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

#include "module_deconv.h"

module_deconv::module_deconv() = default;

std::string module_deconv::get_name()
{
    return "Deconv";
}

void module_deconv::run( options *opts )
{
    options_deconv *d_opts = ( options_deconv * ) opts;

    auto map = parse_linked_file( d_opts->linked_fname );

}

parallel_map<std::size_t, int>
module_deconv::parse_linked_file( std::string fname )
{
    std::ifstream in_stream( fname );

    parallel_map<std::size_t, int>
        ret_map;
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

                            for( auto it = comma_delimited.begin(); it != comma_delimited.end(); ++it )
                                {
                                    std::size_t id = boost::lexical_cast<std::size_t>( *it );

                                    if( ret_map.find( id ) == ret_map.end() ) 
                                        {
                                            ret_map[ id ] = 0;
                                        }
                                    ret_map[ id ] += get_score( comma_delimited.size() );
                                }
                        }
                }
            ++lineno;
        }

    return ret_map;
}

int module_deconv::get_score( std::size_t size )
{
    return 1;
}
