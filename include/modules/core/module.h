#ifndef MODULE_HH_INCLUDED
#define MODULE_HH_INCLUDED

#include "options.h"

class module
{

 public:

    /**
     * Get the name of a module.
     **/
    std::string name;

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

    /**
     * Class to retrieve the name of a module. 
     * @returns string name of a module.
     **/
    virtual std::string get_name();

};
#endif 
