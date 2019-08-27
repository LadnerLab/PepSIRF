#include "module_normalize.h"
#include "options_normalize.h"

module_normalize::module_normalize()
{
    name = "Normalize";
}

std::string module_normalize::get_name()
{
    return name;
}

void module_normalize::run( options *opts )
{
    options_normalize *n_opts = (options_normalize*) opts;

}
