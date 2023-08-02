#ifndef LOGGER_H
#define LOGGER_H

#include <fstream>
#include <iostream>
#include <string>

#include "options.h"

class Log
{
private:
    static std::ofstream logstream;

public:
    Log();

    ~Log();

    static void close();

    static void error(const std::string err_str);

    static void info(const std::string info_str);

    static void open(const std::string logfile);

    static void warn(const std::string warn_str);
};


#endif /* LOGGER_H */

