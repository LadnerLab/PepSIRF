#ifndef OPTIONS_DECONV_HH_INCLUDED
#define OPTIONS_DECONV_HH_INCLUDED

#include <string>
#include "options.h"

class options_deconv : public options
{
 public:
    
    /**
     * Returns the arguments that are stored by 
     * the options object.
     **/
    std::string get_arguments();
};

#endif // OPTIONS_DECONV_HH_INCLUDED
