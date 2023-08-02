#ifndef LOGGER_H
#define LOGGER_H

#include <fstream>
#include <iostream>
#include <string>

#include "options.h"

class Log
{
private:
    static std::ofstream* logstream;

public:
    Log();

    ~Log();

    static void error(std::string err_str);

    static void info(std::string info_str);

    static void warn(std::string warn_str);
};


#endif /* LOGGER_H */

