#include "module_factory.h"

module_factory::module_factory() = default; // default constructor

module *module_factory::create( const char *mod_name )
{
    std::string mod_str = mod_name;

    if( !mod_str.compare( "demux" ) )
        {
            return new module_demux();
        }
    return nullptr;
}
