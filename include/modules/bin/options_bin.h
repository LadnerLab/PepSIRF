#ifndef OPTIONS_BIN_HH_INCLUDED
#define OPTIONS_BIN_HH_INCLUDED

#include "options.h"

class options_bin : public options
{
 public:
    options_bin();
    std::string get_arguments();

};

#endif // OPTIONS_BIN_HH_INCLUDED
