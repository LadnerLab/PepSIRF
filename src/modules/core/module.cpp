#include <iostream>
#include "module.h"

module::module()
{
    name = "default module";
}
module::~module() = default;

void module::run( options *opts )
{
    if( opts->num_threads == 0 )
        {
            std::cout << "Zero threads!\n";
        }
}

std::string module::get_name()
{
    return name;
}
