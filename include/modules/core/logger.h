#ifndef LOGGER_H
#define LOGGER_H

#include <fstream>
#include <iostream>
#include <string>

#include "options.h"

void info(std::string info_str)
{
    std::ofstream log(options::get_logfile());

    std::cout << info_str << "\n";
    log << info_str << "\n";

    log.close();
}


#endif /* LOGGER_H */

