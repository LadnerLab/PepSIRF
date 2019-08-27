#ifndef OPTIONS_NORMALIZE_HH_INCLUDED 
#define OPTIONS_NORMALIZE_HH_INCLUDED 

#include <string>

#include "options.h"

class options_normalize : public options
{
 public:
    options_normalize();

    std::string get_arguments();


};

#endif // OPTIONS_NORMALIZE_HH_INCLUDED
