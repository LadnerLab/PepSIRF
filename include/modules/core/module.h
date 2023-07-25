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
     * Name of the module
     **/
    static std::string name;

    /**
     * Retrieves the name of the module. 
     * @returns string name of the module.
     **/
    static std::string get_name();

    /**
     * Run the module.
     * @param opts Pointer to an options object, each module will have 
     *        its own type of options. 
     **/
    virtual void run( options *opts );
};


#endif 

