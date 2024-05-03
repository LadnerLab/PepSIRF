#ifndef MODULE_INFO_HH_INCLUDED
#define MODULE_INFO_HH_INCLUDED

#include "module.h"
#include "options_info.h"

class module_info : public module
{
public:
    /**
     * Default constructor
     */
    module_info();

    /**
     * Run info module
     * @param options provided by user
     */
    void run( options *opts );
};


#endif /* MODULE_INFO_HH_INCLUDED */

