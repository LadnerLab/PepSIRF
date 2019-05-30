#ifndef MODULE_FACTORY_HH_INCLUDED
#define MODULE_FACTORY_HH_INCLUDED
#include "module.h"
#include "module_demux.h"

/**
 * Class for creating modules based on a string identifier passed as a 
 * parameter to the create method.
 **/
class module_factory
{
 public:
    /**
     * Default constructor
     **/ 
    module_factory();

    /**
     * Method to create a moduel given the name of the module.
     * @param module_name String containing the name of module to create.
     * @pre module_name is a valid module name
     * @note The returned module will be created with 'new', and must be deleted.
     * @returns pointer to a module that has been dynamically allocated, or null if 
     *          an invalid module name is passed
     **/
    module *create( const char *module_name );


};

#endif // MODULE_FACTORY_HH_INCLUDED
