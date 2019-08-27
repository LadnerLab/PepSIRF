#ifndef MODULE_NORMALIZE_HH_INCLUDED
#define MODULE_NORMALIZE_HH_INCLUDED

#include <string>

#include "module.h"
#include "options_normalize.h"


class module_normalize : public module
{

 public:
    std::string name;

    module_normalize();

    std::string get_name();


};

#endif // MODULE_NORMALIZE_HH_INCLUDED
