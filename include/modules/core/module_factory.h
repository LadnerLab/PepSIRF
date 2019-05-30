#ifndef MODULE_FACTORY_HH_INCLUDED
#define MODULE_FACTORY_HH_INCLUDED
#include "module.h"

class module_factory
{
 public:
    module_factory();
    virtual module *create( const char *module_name );


};

#endif // MODULE_FACTORY_HH_INCLUDED
