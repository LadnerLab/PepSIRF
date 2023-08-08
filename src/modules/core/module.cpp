#include <iostream>
#include "module.h"

module::module()
{
    name = "default module";
}

module::~module() = default;

std::string module::name = "";

void module::run( options *opts )
{
    if( opts->num_threads == 0 )
        {
            Log::info("Zero threads!");
        }
}

std::string module::get_name()
{
    return name;
}

