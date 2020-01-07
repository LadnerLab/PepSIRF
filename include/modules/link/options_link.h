#ifndef OPTIONS_LINK_HH_INCLUDED
#define OPTIONS_LINK_HH_INCLUDED
#include "options.h"
#include <string>

class options_link : public options
{
 public:
    options_link();

    std::string get_arguments();

};

#endif // OPTIONS_LINK_HH_INCLUDED
