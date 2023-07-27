#include "logger.h"

void info(std::string info_str)
{
    std::ofstream log(options::get_logfile());

    std::cout << info_str;
    log << info_str;

    log.close();
}

