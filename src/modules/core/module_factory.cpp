#include "module_factory.h"
#include "module_deconv.h"
#include "module_normalize.h"

module_factory::module_factory() = default; // default constructor

module *module_factory::create( const char *mod_name )
{
    std::string mod_str = mod_name;

    if( !mod_str.compare( "demux" ) )
        {
            return new module_demux();
        }
    else if( !mod_str.compare( "deconv" ) )
        {
            return new module_deconv();
        }
    else if( !mod_str.compare( "normalize" ) )
        {
            return new module_normalize();
        }


    return nullptr;
}
