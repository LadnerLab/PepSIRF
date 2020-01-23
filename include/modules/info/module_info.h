#ifndef MODULE_INFO_HH_INCLUDED
#define MODULE_INFO_HH_INCLUDED
#include "module.h"
#include "options_info.h"

class module_info : public module
{
public:
    module_info() = default;
    void run( options *opts );

};

#endif // MODULE_INFO_HH_INCLUDED
