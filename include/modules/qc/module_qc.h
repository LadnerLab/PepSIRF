#include "fif_parser.h"
#include "module.h"
#include "options.h"

class module_qc : module  
{
    public:
        std::string name;
        module_qc();
        void run ( options *opts );
        std::string get_name( );
};