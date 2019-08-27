#include <algorithm>
#include "options_factory.h"
#include "options_deconv.h"
#include "options_normalize.h"


options *options_factory::create( int argc, char ***argv )
{
    std::string mod_name;
    // make sure at least a module name was supplied
    if( argc >= 2 )
        {
            mod_name = (*argv)[ 1 ];

            // make sure that the module name is lower-case
            std::transform( mod_name.begin(), mod_name.end(), mod_name.begin(), ::tolower );

            if( !mod_name.compare( "demux" ) )
                {
                    return new options_demux();
                }
            else if( !mod_name.compare( "deconv" ) )
                {
                    return new options_deconv();
                }
            else if( !mod_name.compare( "normalize" ) )
                {
                    return new options_normalize();
                }
        }

    return nullptr;
}
