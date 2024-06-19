#include "module.h"
#include "options.h"

options::~options() = default;

std::string options::logfile = "";

std::string options::get_arguments()
{
    return "Module name";
}

std::string options::set_default_log()
{
    std::time_t exec_time = std::time(nullptr);
    std::string local_time = std::asctime(std::localtime(&exec_time));
    local_time = local_time.substr(11, 8);

    size_t pos = 0;
    while ((pos = local_time.find(":", pos)) != std::string::npos)
    {
        local_time.replace(pos, 1, "-");
    }

    return module::get_name() + "_" + local_time + ".log";
}

std::string options::get_logfile()
{
    return logfile;
}

