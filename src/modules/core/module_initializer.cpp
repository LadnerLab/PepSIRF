#include "module_initializer.h"
#include <algorithm>

module_initializer::module_initializer() 
    : opts( nullptr ), mod( nullptr ), opt_parser( nullptr )
      {} 

module_initializer::~module_initializer() 
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


void module_initializer::initialize( const std::string& mod_name )
{

    if( !mod_name.compare( "subjoin" ) )
        {
            opts        = new options_subjoin();
            opt_parser  = new options_parser_subjoin();
            mod         = new module_subjoin();
        }
    else
        {
            throw std::runtime_error( "Invalid module name entered." );
        }
}

options *module_initializer::get_opts()
{
    throw_on_null( opts );

    return opts;
}

module *module_initializer::get_module()
{
    throw_on_null( mod );

    return mod;
}

options_parser *module_initializer::get_options_parser()
{
    throw_on_null( opt_parser );

    return opt_parser;
}

