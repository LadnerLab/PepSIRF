#ifndef MODULE_P_ENRICH_HH_INCLUDED
#define MODULE_P_ENRICH_HH_INCLUDED
#include "module.h"
#include "options_p_enrich.h"

class module_p_enrich : public module
{

public:

    /**
     * Run the bin module, passing 
     * the options specified by a user.
     **/
    void run( options *opts );


};

#endif // MODULE_P_ENRICH_HH_ENCLUDED
