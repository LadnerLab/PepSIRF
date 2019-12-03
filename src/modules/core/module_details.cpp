#include "module_details.h"

module_details::module_details() 
    : opts( nullptr ), mod( nullptr ), opt_parser( nullptr ),
      initialized( false ) {} 

module_details::~module_details() 
{
    if( initialized )
        {
            delete opts;
            delete mod;
            delete opt_parser;
        }
}


void module_details::initialize( const std::string& module_name )
{

    initialized = true;
}
