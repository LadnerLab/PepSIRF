#include "options.h"

options::~options() = default;
std::string options::logfile = "";

std::string options::get_arguments()
{
    return "Module name";
}

std::string options::get_default_log()
{
    return  "module_name.log";
}

std::string options::get_logfile()
{
    return logfile;
}


