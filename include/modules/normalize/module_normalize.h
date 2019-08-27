#ifndef MODULE_NORMALIZE_HH_INCLUDED
#define MODULE_NORMALIZE_HH_INCLUDED

#include <string>

#include "module.h"
#include "options.h"

/**
 * Class for the PepSIRF normalization module.
 **/
class module_normalize : public module
{

 public:
    std::string name;

    module_normalize();

    std::string get_name();

    void run( options *opts );

};

#endif // MODULE_NORMALIZE_HH_INCLUDED
