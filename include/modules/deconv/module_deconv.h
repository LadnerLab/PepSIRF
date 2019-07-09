#ifndef MODULE_DECONV_HH_INCLUDED 
#define MODULE_DECONV_HH_INCLUDED 
#include <string>

#include "module.h"
#include "options_deconv.h"

class module_deconv : public module
{
 public:
    module_deconv();
    std::string get_name();
    void run( options *opts );
};

#endif // MODULE_DECONV_HH_INCLUDED 
