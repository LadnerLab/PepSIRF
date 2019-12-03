#include "module_details.h"
#include <algorithm>

module_details::module_details() 
    : opts( nullptr ), mod( nullptr ), opt_parser( nullptr )
      {} 

module_details::~module_details() 
{
    if( opts )
        {
            delete opts;
        }
    if( mod )
        {
            delete mod;
        }
    if( opt_parser )
        {
            delete opt_parser;
        }
}


void module_details::initialize( const std::string& mod_name )
{

    if( !mod_name.compare( "demux" ) )
        {
            opts       = new options_demux();
            opt_parser = new options_parser_demux();
            mod        = new module_demux();
        }
    else if( !mod_name.compare( "deconv" ) )
        {
            opts       = new options_deconv();
            opt_parser = new options_parser_deconv();
            mod        = new module_deconv();
        }
    else if( !mod_name.compare( "norm" ) )
        {
            opts        = new options_normalize();
            opt_parser  = new options_parser_normalize();
            mod         = new module_normalize();
        }
    else
        {
            throw std::runtime_error( "Invalid module name entered." );
        }
}

options *module_details::get_opts()
{
    throw_on_null( opts );

    return opts;
}

module *module_details::get_module()
{
    throw_on_null( mod );

    return mod;
}

options_parser *module_details::get_options_parser()
{
    throw_on_null( opt_parser );

    return opt_parser;
}

