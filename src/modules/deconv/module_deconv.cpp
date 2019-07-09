#include "module_deconv.h"
#include <iostream>

module_deconv::module_deconv() = default;

std::string module_deconv::get_name()
{
    return "Deconv";
}

void module_deconv::run( options *opts )
{
    std::cout << "Module deconv!\n";
}
