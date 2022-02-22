#include "fif_parser.h"
#include "module.h"
#include "options.h"

class module_qc : public module  
{
    public:
        std::string name;
        module_qc() = default;
        void run ( options *opts );
        std::string get_name( );
        std::vector<std::string> get_combined( std::size_t iteration, std::vector<std::vector<std::string>> indata, std::vector<std::string> outdata);
};