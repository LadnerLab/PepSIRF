#ifndef OPTIONS_SUBJOIN_HH_INCLUDED
#define OPTIONS_SUBJOIN_HH_INCLUDED
#include "options.h"

/**
 * Options for the subjoin module of the 
 * PepSIRF package.
 **/
class options_subjoin : public options
{
    std::string get_arguments();
};

#endif // OPTIONS_SUBJOIN_HH_INCLUDED
