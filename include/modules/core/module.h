#ifndef MODULE_HH_INCLUDED
#define MODULE_HH_INCLUDED

#include "options.h"

class module
{

 public:
    /**
     * Default constructor.
     **/
    module(); 

    /**
     * Default deconstructor
     **/
    virtual ~module();

    /**
     * Run a module.
     * @param opts Pointer to an options object, each module will have 
     *        its own type of options. 
     **/
    virtual void run( options *opts );

};
#endif 
