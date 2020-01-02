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
    else if( !mod_name.compare( "subjoin" ) )
        {
            opts        = new options_subjoin();
            opt_parser  = new options_parser_subjoin();
            mod         = new module_subjoin();
        }
    else if( !mod_name.compare( "zscore" ) )
        {
            opts        = new options_zscore();
            opt_parser  = new options_parser_zscore();
            mod         = new module_zscore();
        }
    else if( !mod_name.compare( "bin" ) )
        {
            opts        = new options_bin();
            opt_parser  = new options_parser_bin();
            mod         = new module_bin();
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

