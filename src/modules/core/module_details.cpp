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


void module_details::initialize( const std::string& module_name )
{

    initialized = true;
}
