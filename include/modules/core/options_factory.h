#ifndef OPTIONS_FACTORY_HH_INCLUDED 
#define OPTIONS_FACTORY_HH_INCLUDED 
#include "options.h"
#include "options_factory.h"

class options_factory
{
 public:
    options *create( int argc, char ***argv );
};

#endif // OPTIONS_FACTORY_HH_INCLUDED 
