#include <iostream>
#include "module.h"

module::module() = default;
module::~module() = default;

void module::run( options *opts )
{
    if( opts->num_threads == 0 )
        {
            std::cout << "Zero threads!\n";
        }
}
