#include <vector>
#include "modules.h"
#include <iostream>

std::string modules::join_names( const std::string& on )
{
    return join_names_from_iterable( module_names, on );
}


std::string modules::join_names( const std::string& prefix,
                                 const std::string& on
                               )
{
    std::vector<std::string> mod_names;

    for( const auto& mod_name : module_names )
        {
            std::string current( "" );

            current = prefix + mod_name;
            mod_names.emplace_back( current );
        }

    return modules::join_names_from_iterable( mod_names, on );
}


