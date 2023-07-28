#include "logger.h"

// TODO: deal with repeated opening and closing log file
void error(std::string err_str)
{
    std::ofstream log(options::get_logfile());

    std::cout << "Error: " << err_str;
    log << "Error: " << err_str;

    log.close();
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
    std::ofstream log(options::get_logfile());

    std::cout << "Warning: " << warn_str;
    log << "Warning: " << warn_str;

    log.close();
}

