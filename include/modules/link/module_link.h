#ifndef MODULE_LINK_HH_INCLUDED 
#define MODULE_LINK_HH_INCLUDED 
#include "module.h"
#include "options_link.h"

class module_link : public module
{
 public:
    module_link();

    /**
     * Run the link module, passing the options
     * specified by a user.
     **/
    void run( options *opts );

};

#endif // MODULE_LINK_HH_INCLUDED
