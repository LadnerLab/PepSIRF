#include "logger.h"

void error(std::string err_str)
{
}

void info(std::string info_str)
{
    std::ofstream log(options::get_logfile());

    std::cout << info_str;
    log << info_str;

    log.close();
}

void warn(std::string warn_str)
{
}

