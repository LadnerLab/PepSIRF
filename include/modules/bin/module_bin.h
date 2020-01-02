#ifndef MODULE_BIN_HH_INCLUDED
#define MODULE_BIN_HH_INCLUDED

#include "module.h"
#include "options_bin.h"

class module_bin : public module
{
 public:
    module_bin();

    /**
     * Run the bin module, passing 
     * the options specified by a user.
     **/
    void run( options *opts );


};

#endif // MOUDLE_BIN_HH_INCLUDED
