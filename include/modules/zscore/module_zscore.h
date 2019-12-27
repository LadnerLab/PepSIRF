#ifndef MODULE_ZSCORE_HH_INCLUDED
#define MODULE_ZSCORE_HH_INCLUDED
#include "module.h"
#include "options_zscore.h"

class module_zscore : public module
{
 public:
    module_zscore();
    void run( options *opts );

};

#endif // MODULE_ZSCORE_HH_INCLUDED
