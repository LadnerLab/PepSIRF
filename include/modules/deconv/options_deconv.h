#ifndef OPTIONS_DECONV_HH_INCLUDED
#define OPTIONS_DECONV_HH_INCLUDED

#include "options.h"

#include <string>
#include <vector>

class options_deconv : public options
{
 public:
    
    /**
     * Returns the arguments that are stored by 
     * the options object.
     **/
    std::string get_arguments();

    std::string linked_fname;
    std::size_t threshold;
};

#endif // OPTIONS_DECONV_HH_INCLUDED
