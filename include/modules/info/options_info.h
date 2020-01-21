#ifndef OPTIONS_INFO_HH_INCLUDED
#define OPTIONS_INFO_HH_INCLUDED

#include <string>
#include "options.h"

class options_info : public options
{

public:
    options_info() = default;

    std::string get_arguments();

};

#endif // OPTIONS_INFO_HH_INCLUDED
