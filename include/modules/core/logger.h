#ifndef LOGGER_H
#define LOGGER_H

#include <iostream>
#include <fstream>
#include "options.h"

using namespace std;

void error(const char* format)
{
}

void warn(const char* format)
{
}

void info(const char* format)
{
    ofstream log(options::get_logfile());

    std::cout << format << "\n";
    log << format << "\n";

    log.close();
}

#endif /* LOGGER_H */

